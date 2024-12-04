//Replace each extension at each process
def replaceExtension(path, newExtension) {
    return path.getSimpleName().split("_")[0] + newExtension
    //return path.getSimpleName() + newExtension
}

//Convert the demultiplexed, raw sequencing FASTQ files to BAM
//
process convertFastqToSam {
    // publishDir allow you to define the path of output
    publishDir params.outdir

    //the input received as a pairs of FASTQ files
    //
    input:
    tuple val('sample_id'), path(fastq)

    output:
    path "${replaceExtension(fastq[0], '_unmapped.bam')}"

    """
    gatk FastqToSam \
    --FASTQ ${fastq[0]} \
    --FASTQ2 ${fastq[1]} \
    --OUTPUT "${replaceExtension(fastq[0], '_unmapped.bam')}" \
    --SAMPLE_NAME Mysample1
    """
}

// 1.  extraction of UMIs from the insert reads
// it is the parameter "-r" that define the number of intial bases to extract for UMIs. in this pipeline we extract the initial 3 bases for UMIs
// 5M2S+T 5M2S+T : in case of Twist kit
//
process ExtractUmis {
    tag "${sample}: ${convertFastqToSamOut.first()}"
    publishDir params.outdir

    input:
    path convertFastqToSamOut

    output:
    path "${replaceExtension(convertFastqToSamOut, '_umi_extracted.bam')}"

    script:
    def bam = convertFastqToSamOut.first()


    """
    fgbio ExtractUmisFromBam \\
    -i ${convertFastqToSamOut} \
    -o ${replaceExtension(convertFastqToSamOut, '_umi_extracted.bam')} \
    -r 5M2S+T 5M2S+T \
    -t RX \
    -a true
    """
}

// Convert the BAM file with UMI extracted reads to a FASTQ file
//
process convertSamToFastq {
    publishDir params.outdir

    input:
    path ExtractUmisOut

    //two output fastq file
    output:
    tuple val('sample_id'), path("*.fastq")

    // "".baseName.take(30)" allow to take only the 30 first lettre of the file name.
    """
    gatk SamToFastq \
    -I ${ExtractUmisOut} \
    -F "${replaceExtension(ExtractUmisOut, '_R1.fastq')}" \
    -F2 "${replaceExtension(ExtractUmisOut, '_R2.fastq')}" \
    --CLIPPING_ATTRIBUTE XT \
    --MAX_RECORDS_IN_RAM 50000000 \
    --CLIPPING_ACTION 2

    """
}

// adapter and quality trimming
//
process fastpp {
    publishDir params.outdir

    input:
    tuple val('sample_id'), path(convertSamToFastqOut)


    output:
    tuple val('sample_id'), path("*.fastq*")

    script:
    """
    fastp \\
    -i ${convertSamToFastqOut[0]} \\
    -o "${replaceExtension(convertSamToFastqOut[0], '_umi_extracted_trimmed_R1.fastq')}" \\
    -I ${convertSamToFastqOut[1]} \\
    -O "${replaceExtension(convertSamToFastqOut[1], '_umi_extracted_trimmed_R2.fastq')}" \\
    -g -W 5 -q 20 -u 40 -x -3 -l 75 -c \
    -j fastp.json \
    -h fastp.html \
    -w 12
    """
}

//  Alignement before deduplication: Align quality and adapter trimmed reads to the reference genome
//
process bwaMem {
    tag "${sample}: ${fastppOut.first()}"
    publishDir params.outdir

    input:
    path params.ref
    tuple val(sample), path(fastppOut)

    output:
    file "${replaceExtension(fastppOut[0], '_umi1_extracted_aligned.bam')}"

    script:
    def bam = fastppOut.first()

    """
    bwa mem \
    -t 10 \
    -M ${params.ref} \\
    ${fastppOut[0]} \\
    ${fastppOut[1]} \\
    | \
    samtools view -Sb -o ${replaceExtension(fastppOut[0], '_umi1_extracted_aligned.bam')}
    """
}

// Extraction of read quality metrics before deduplication
//
process collectmetrics1 {
    publishDir params.outdir

    input:
    path params.ref // reference genome
    path interval_list //
    path bwaMemOut  // the output of bwaMem process

    output:
    file "${replaceExtension(bwaMemOut, '_output_hs_metrics1.txt')}"


    """
    picard CollectHsMetrics \
    --REFERENCE_SEQUENCE $params.ref \
    --BAIT_INTERVALS ${interval_list} \
    --TARGET_INTERVALS ${interval_list} \
    --INPUT ${bwaMemOut} \
    --OUTPUT "${replaceExtension(bwaMemOut, '_output_hs_metrics1.txt')}"
    """
}

// Generate a multi-quality control report from collected metrics data (collectmetrics1 output).

 process multiQc1{
     publishDir params.outdir

     input:
     path collectmetrics1Out  // the output of collectmetrics1 process

     output:
     file "${replaceExtension(collectmetrics1Out, '_reportebefore.html')}"

     """
     multiqc ${collectmetrics1Out} -o "${replaceExtension(collectmetrics1Out, '_reportebefore.html')}"
     """
 }

// Merge the two BAM files containing:
//1: the UMI information: output of ExtractUmis process
//2: the alignment coordinate information: output of bwaMEM process
//
process MergeBam {
    tag "${sample}: ${umi_extracted_aligned[0].first()}: ${umi_extracted_aligned[1].name}"
    publishDir params.outdir

    // The outcomes from both the ExtractUmis process and the bwaMEM process are merged according to their filenames into lists for each sample
    input:
    tuple val(sample), path(umi_extracted_aligned)
    path params.ref

    output:
    file "${replaceExtension(umi_extracted_aligned[0], '_umi_extracted_aligned_merged.bam')}"

    script:
    def bam = umi_extracted_aligned[0].first()

    """
    gatk MergeBamAlignment \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_REMOVE NM \
    --ATTRIBUTES_TO_REMOVE MD \
    --ALIGNED_BAM ${umi_extracted_aligned[0]} \
    --UNMAPPED_BAM ${umi_extracted_aligned[1]} \
    --OUTPUT "${replaceExtension(umi_extracted_aligned[0], '_umi_extracted_aligned_merged.bam')}" \
    --REFERENCE_SEQUENCE ${params.ref} \
    --SORT_ORDER 'queryname' \
    --ALIGNED_READS_ONLY true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --CLIP_OVERLAPPING_READS false
    """
}



process umiMergeFilt {
    publishDir params.outdir

    input:
    path MergeBamOut

    output:
    file "${replaceExtension(MergeBamOut, '_merged_filtered.bam')}"

    script:
    """
    samtools view \
    -f2 \
    -bh ${MergeBamOut} \
    > \
    ${replaceExtension(MergeBamOut, '_merged_filtered.bam')}
    """
}

// 2. Identify and group reads originating from the same source molecule
// The user can control how many errors/mismatches are allowed in the UMI sequence when assigning source molecules (--edits=n).
//
process GroupReads {
    publishDir params.outdir

    input:
    path umiMergeFiltOut


    output:
    file "${replaceExtension(umiMergeFiltOut, '_umi_grouped.bam')}"

    script:
    """
    fgbio GroupReadsByUmi \\
    -i ${umiMergeFiltOut} \\
    -o '${replaceExtension(umiMergeFiltOut, '_umi_grouped.bam')}' \\
    --strategy=adjacency \
    --edits=1 \
    -t RX \
    -f family_size_counts.txt
    """
}

// 3. Generate consensus reads
// Calculate the consensus sequence, Reads that occur as singletons are discarded by default but this can be changed
// by setting the –min-reads flag to 1, in so doing the single read will be considered the consensus.
//
process CallConsensus {
    tag "${sample}: ${GroupReadsOut.first()}"

    publishDir params.outdir

    input:
    path GroupReadsOut


    output:
    file "${replaceExtension(GroupReadsOut, '_consensus_unmapped.bam')}"

    script:
    def bam = GroupReadsOut.first()

    """
    fgbio CallMolecularConsensusReads \\
    -i ${GroupReadsOut} \\
    -o '${replaceExtension(GroupReadsOut, '_consensus_unmapped.bam')}' \\
    --error-rate-post-umi 40 \
    --error-rate-pre-umi 45 \
    --output-per-base-tags false \
    --min-reads 2 \
    --max-reads 50 \
    --read-name-prefix='consensus'
    """
}

// Following consensus calling, the UMI group collapsing process leads to the removal of alignment coordinate information.
// To address this, it is necessary to convert the consensus_unmapped.bam file to the FASTQ format.
//
process covertSamToFastq {
    publishDir params.outdir

    input:
    path CallConsensusOut

    //two output fastq file
    output:
    tuple val('sample_id'), path("*.fastq*")

    """
    gatk SamToFastq \
    -I ${CallConsensusOut} \
    --F "${replaceExtension(CallConsensusOut, '_consensus_unmapped_R1.fastq')}" \
    --F2 "${replaceExtension(CallConsensusOut, '_consensus_unmapped_R2.fastq')}" \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2 \

    """
}

// Align the FASTQ files to the reference with BWA-mem and output to a new consensus_mapped.bam file
//
process bwaMem2 {
    tag "${sample}: ${covertSamToFastqOut.first()}"
    publishDir params.outdir

    input:
    path params.ref
    tuple val(sample), path(covertSamToFastqOut)

    output:
    file "${replaceExtension(covertSamToFastqOut[0], '_consensus_mapped.bam')}"

    script:
    def bam = covertSamToFastqOut.first()

    """
    bwa mem \
    -v 3 \
    -t 10 \
    -Y \
    -M \
    ${params.ref} \\
    ${covertSamToFastqOut[0]} \\
    ${covertSamToFastqOut[1]} \\
    | \
    samtools view -bh - > ${replaceExtension(covertSamToFastqOut[0], '_consensus_mapped.bam')}
    """
}

// Extraction of read quality metrics after deduplication
//
process collectmetrics2 {
    publishDir params.outdir

    input:
    path params.ref
    path interval_list
    path bwaMem2Out

    output:
    file "${replaceExtension(bwaMem2Out, '_output_hs_metrics2.txt')}"


    """
    picard CollectHsMetrics \
    --REFERENCE_SEQUENCE $params.ref \
    --BAIT_INTERVALS ${interval_list} \
    --TARGET_INTERVALS ${interval_list} \
    --COVERAGE_CAP 1000 \
    --INPUT ${bwaMem2Out} \
    --OUTPUT ${replaceExtension(bwaMem2Out, '_output_hs_metrics2.txt')}
    """
}

// Generate a multi-quality control report from collected metrics data (process collectmetrics2 output).
//
process multiQc{
    publishDir params.outdir
    input:
    path collectmetrics2Out

    output:
    file "${replaceExtension(collectmetrics2Out, '_reportAfter')}"

    """
    multiqc --force ${collectmetrics2Out} -o "${replaceExtension(collectmetrics2Out, '_reportAfter')}"
    """
}

// Sort the consensus_mapped.bam with the consensus_unmapped.bam to prepare them as input for the next step
//
process sortConsensus2 {
    tag "${sample}: ${consensus_mapped1[0].first()}: ${consensus_mapped1[1].name}"
    publishDir params.outdir

    input:
    tuple val(sample), path(consensus_mapped1)

    output:
    tuple val('sample'), path("*.bam")

    script:

     """
    gatk SortSam \
    -I ${consensus_mapped1[0]} \
    --OUTPUT "${replaceExtension(consensus_mapped1[0], '_consensus_mapped_sorted.bam')}" \
    --SORT_ORDER queryname

    gatk SortSam \
    -I ${consensus_mapped1[1]} \
    --OUTPUT "${replaceExtension(consensus_mapped1[1], '_consensus_unmapped_sorted.bam')}" \
    --SORT_ORDER queryname
    """
}

// Finally, merge the consensus_mapped.bam with the consensus_unmapped.bam to retain the UMI group information.
//
process MergeBam2 {
    tag "${sample}: ${sortConsensus2Out[0].first()}"

    publishDir params.outdir

    input:
    tuple val('sample_id'), path(sortConsensus2Out)
    path params.ref

    output:
    file "${replaceExtension(sortConsensus2Out[0], '_consensusMerge.bam')}"

    """
    gatk MergeBamAlignment \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_RETAIN RX \
    --ALIGNED_BAM ${sortConsensus2Out[0]} \
    --UNMAPPED_BAM ${sortConsensus2Out[1]} \
    --OUTPUT "${replaceExtension(sortConsensus2Out[0], '_consensusMerge.bam')}" \
    --REFERENCE_SEQUENCE ${params.ref} \
    --SORT_ORDER coordinate \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --CLIP_OVERLAPPING_READS false
    """
}

// Generate the index (.bai) of the MergeBam2 output bam file
//
process bamindex {
    tag "${sample}: ${MergeBam2Out.first()}"

    publishDir params.outdir

    input:
    path MergeBam2Out

    output:
    path "${replaceExtension(MergeBam2Out, '_consensusMerge.bam.bai')}"

    """
    gatk BuildBamIndex \
    --INPUT ${MergeBam2Out} \
    --OUTPUT "${replaceExtension(MergeBam2Out, '_consensusMerge.bam.bai')}"
    """
}

// Variant calling step using vardict
//
process varCallvardict {
    tag "${sample}: ${bam[0].first()}: ${bam[1].name}"

    publishDir params.outdir

   input:
    path params.ref
    tuple val('sample_id'), path(bam)
    path bed

    output:
    file "${replaceExtension(bam[0], '_consensus_vardict.vcf')}"

    """
    vardict -G $params.ref -f 0.0005 -N sample_name -b ${bam[0]} -c 1 -S 2 -E 3 -g 4 $bed | $params.teststrandbias | $params.var2vcf > ${replaceExtension(bam[0], '_consensus_vardict.vcf')}
    """
}

// Annotation step using VEP
//
process annotationVep{
    publishDir params.outdir

    input:
    path varCallvardictOut
    path cach
    path fasta

    output:
    file "${replaceExtension(varCallvardictOut, '_vep.vcf')}"

    """
    vep -i ${varCallvardictOut} -o ${replaceExtension(varCallvardictOut, '_vep.vcf')} --cache --dir_cache ${cach} --offline --force_overwrite --vcf --numbers --refseq --symbol --hgvs --canonical --max_af --fasta ${fasta}

    """
}

// Transforming VCF file of annotation step into tsv format
//
process vcf2tsv{
    publishDir params.outdir
    input:
    path annotationVepOut

    output:
    file "${replaceExtension(annotationVepOut, '.tsv')}"

    """
    vcf2tsvpy --input_vcf ${annotationVepOut} --out_tsv ${replaceExtension(annotationVepOut, '.tsv')}
    """
}


workflow {

// step 1: Extraction of UMIs from the insert reads

    // in this step the input received as a pairs of FASTQ files so we transfomr each pair to a tuple  named read_pairs_fastq
   Channel
            .fromFilePairs(params.fastq,checkIfExists:true)
            .set {read_pairs_fastq}
    convert1 = convertFastqToSam(read_pairs_fastq)
    extract1 = ExtractUmis (convert1)
    convert2 = convertSamToFastq(extract1)
    FASTP= fastpp(convert2)
    bwa1= bwaMem(params.ref, FASTP)
    collectMet1 = collectmetrics1(params.ref, params.bedInterval, bwa1)
    // multiQc1(collectMet1)

// Combine the two channels (extract1 and bwa1 ) and group by file name

    bwaMem.out
        .set { bwa1_out } // name the bwaMem output as bwa1_out
    ExtractUmis.out
        .set { extract1_out } // name the ExtractUmis output as extract1_out

    bwa1_out
        .map { filepath -> [filepath.name.tokenize('_')[0], filepath] }
        // Extract the base sample name (e.g., 'P22D0165ct') and pair it with the file path
        .set { bwa1_out_tuple }

    extract1_out
        .map { filepath -> [filepath.name.tokenize('_')[0], filepath] }
        // Extract the base sample name (e.g., 'P22D0165ct') and pair it with the file path
        .set { extract1_out_tuple }

    // Combine and group bwa1_out_tuple and extract1_out_tuple based on the shared sample name
    bwa1_out_tuple
        .mix(extract1_out_tuple) // Combine the two channels
        .groupTuple() // Group by the first element of the tuple (sample name)
        .flatMap { sample, path_list ->
            path_list.split {
                it.name.endsWith('extracted_aligned.bam') // Split by files ending with 'extracted_aligned.bam'
            }.combinations()
        }
        .map { combination ->
            def fmeta = [ "id": "test" ] // Add metadata (adjust as needed)
            [ fmeta, combination*.first(), combination*.last() ] // Include metadata and file paths in output
            [ fmeta.id, combination ]
        }
        .set { input_merge } // Store processed combinations as input_merge

    MergeBam(input_merge, params.ref) // Call the MergeBam process with input_merge

    umiFiltr= umiMergeFilt(MergeBam.out)

// step 2: Group reads by UMI

    groupR= GroupReads(umiFiltr)

// step 3: Generate consensus reads

    consensusR1= CallConsensus(groupR)
    samToFastq= covertSamToFastq(CallConsensus.out)
    bwa2= bwaMem2(params.ref, samToFastq)
    collectMet2 = collectmetrics2(params.ref, params.bedInterval, bwa2)
    multiQc(collectMet2)

    // Combine the two channels (bwa2 and consensusR1 ) and group by file name
    bwaMem2.out
        .set { bwa2_out }
    CallConsensus.out
        .set { consensusR1_out }

    bwa2_out
        .map { filepath -> [filepath.name.tokenize('_')[0], filepath ] }
        .set { bwa2_out_tuple }
    consensusR1_out
        .map { filepath -> [filepath.name.tokenize('_')[0], filepath ] }
        .set { consensusR1_out_tuple }

     bwa2_out_tuple
        .mix(consensusR1_out_tuple)
        .groupTuple()
        .flatMap { sample, path_list ->
            path_list.split {
            it.name.endsWith('_consensus_mapped.bam')
      }.combinations()
    }
    .map { combination ->
        def fmeta = [ "id": "test" ]

        [ fmeta, combination*.first(), combination*.last() ]
        [ fmeta.id, combination ]
    }
    .set { input_sortConsensus2 } // stores these processed combinations as input_sortConsensus2

    sortConsensus2(input_sortConsensus2)

    mergeBam2= MergeBam2(sortConsensus2.out, params.ref)
    indexBam= bamindex(mergeBam2)

    // Combine the two channels (mergeBam2 and indexBam ) and group by file name
    MergeBam2.out
        . set { mergeBam2_out }
    bamindex.out
        . set { indexBam_out }

    mergeBam2_out
        .map { filepath -> [filepath.name.toString().tokenize('.')[0], filepath ] }
        .set { mergeBam2_out_tuple }
    indexBam_out
        .map { filepath -> [filepath.name.toString().tokenize('.')[0], filepath ] }
        .set { indexBam_out_tuple }

     mergeBam2_out_tuple
        .mix(indexBam_out_tuple)
        .groupTuple()
        .flatMap { sample, path_list ->
        path_list.split {
        it.name.endsWith('_consensusMerge.bam')
      }.combinations()
    }
    .map {
        def fmeta = [ "id": "test" ]

        [ fmeta, it*.first(), it*.last() ]
        [ fmeta.id, it ]
    }
    .set { input_vardict } // stores these processed combinations as input_vardict

// variant calling step
    varCallvardict(params.ref, input_vardict, params.bed)
// annotation step
    vepAnn = annotationVep(varCallvardict.out, params.cach, params.fasta)

}
