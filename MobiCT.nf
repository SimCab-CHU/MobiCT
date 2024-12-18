// PROCESSES

// Convert the demultiplexed, raw sequencing FASTQ files to BAM
process ConvertFastqToSam {
    tag "$sample_id"

    input:
        tuple val('sample_id'), path(fastq)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam")

    """
    gatk FastqToSam \
        --FASTQ ${fastq[0]} \
        --FASTQ2 ${fastq[1]} \
        --OUTPUT ${sample_id}${extension}.bam \
        --SAMPLE_NAME ${sample_id} \
        --TMP_DIR ${params.tmp_dir}
    """
}

// Extraction of UMIs from the insert reads
// It is the parameter "-r" that define the number of intial bases to extract
// for UMIs. in this pipeline we extract the initial 3 bases for UMIs
// Example : 5M2S+T 5M2S+T (in case of Twist kit)
process ExtractUmis {
    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam")

    """
    fgbio ExtractUmisFromBam \
        -i ${bam_file} \
        -o ${sample_id}${extension}.bam \
        -r 5M2S+T 5M2S+T \
        -t RX \
        -a true
    """
}

// Convert the BAM file with UMI extracted reads to a FASTQ file
process ConvertSamToFastq {
    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam_file)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.R[1,2].fq")

    """
    gatk SamToFastq \
        -I ${bam_file} \
        -F ${sample_id}${extension}.R1.fq \
        -F2 ${sample_id}${extension}.R2.fq \
        --CLIPPING_ATTRIBUTE XT \
        --CLIPPING_ACTION 2 \
        --MAX_RECORDS_IN_RAM 50000000 
    """
}

// Adapter and quality trimming
process Fastp {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.{json,html}'

    input:
        tuple val(sample_id), path(fastq)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.R[1,2].fq")
        file "${sample_id}${extension}_fastp.json"
        file "${sample_id}${extension}_fastp.html"

    script:
    """
    fastp \
        -i ${fastq[0]} \
        -o ${sample_id}${extension}.R1.fq \
        -I ${fastq[1]} \
        -O ${sample_id}${extension}.R2.fq \
        -g -W 5 -q 20 -u 40 -x -3 -l 75 -c \
        -j ${sample_id}.QC.fastp.json \
        -h ${sample_id}.QC.fastp.html \
        -w 12
    """
}

// Alignement before deduplication: Align quality and adapter trimmed reads to
// the reference genome
process BWAmem {
    tag "$sample_id"
    clusterOptions '-n 10'
    
    input:
        tuple val(sample_id), path(fastq)
        val opt_bwa
        val extension
    
    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam")

    """
    bwa mem \
        ${opt_bwa} \
        -M ${params.ref} \
        ${fastq[0]} \
        ${fastq[1]} \
        | \
        samtools view -bh -o ${sample_id}${extension}.bam
    """
}

// Merge the two BAM files containing:
// 1: the UMI information: output of ExtractUmis process
// 2: the alignment coordinate information: output of bwaMEM process
process MergeBam {
    tag "$sample_id"
    
    input:
        tuple val(sample_id), path(bam_aligned), path(bam_unmapped)
        val extension


    output:
        tuple val(sample_id), file("${sample_id}${extension}.ba[i,m]")

    """
    gatk MergeBamAlignment \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_REMOVE NM \
        --ATTRIBUTES_TO_REMOVE MD \
        --ALIGNED_BAM ${bam_aligned} \
        --UNMAPPED_BAM ${bam_unmapped} \
        --OUTPUT ${sample_id}${extension}.bam \
        --REFERENCE_SEQUENCE ${params.ref} \
        --SORT_ORDER 'queryname' \
        --ALIGNED_READS_ONLY true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --CREATE_INDEX true \
        --CLIP_OVERLAPPING_READS false
    """
}

// 
process UmiMergeFilt {
    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam")

    """
    samtools view \
        -f2 \
        -bh ${bam} \
        > ${sample_id}${extension}.bam
    """
}

// Identify and group reads originating from the same source molecule
// The user can control how many errors/mismatches are allowed in the UMI sequence when assigning source molecules (--edits=n).
process GroupReads {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.txt'

    input:
        tuple val(sample_id), path(bam)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam"), emit: nextout
        file "${sample_id}.family_size_counts.txt"

    """
    fgbio GroupReadsByUmi \
        -i ${bam} \
        -o ${sample_id}${extension}.bam \
        --strategy=adjacency \
        --edits=1 \
        -t RX \
        -f ${sample_id}.QC.family_size_counts.txt
    """
}

// Generate consensus reads
// Calculate the consensus sequence, Reads that occur as singletons are
// discarded by default but this can be changed by setting the –min-reads flag
// to 1, in so doing the single read will be considered the consensus.
process CallConsensus {
    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam")

    """
    fgbio CallMolecularConsensusReads \
        -i ${bam} \
        -o ${sample_id}${extension}.bam \
        --error-rate-post-umi 40 \
        --error-rate-pre-umi 45 \
        --output-per-base-tags false \
        --min-reads 2 \
        --max-reads 50 \
        --read-name-prefix='consensus'
    """
}

workflow RerunConvertSamToFastq {
    take:
        tuple_input
        extension
    main:
        ConvertSamToFastq(tuple_input, extension)
    emit:
        final_out = ConvertSamToFastq.out
}

workflow RerunBWAmem {
    take:
        tuple_input
        opt_bwa
        extension
    main:
        BWAmem(tuple_input, opt_bwa, extension)
    emit:
        final_out = BWAmem.out
}

// Sort the consensus_mapped.bam with the consensus_unmapped.bam to prepare them as input for the next step
process SortConsensus {
    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam)
        val extension

    output:
        tuple val(sample_id), path("${sample_id}${extension}.sort.bam")

    script:

    """
    gatk SortSam \
        -I ${bam} \
        --OUTPUT ${sample_id}${extension}.sort.bam \
        --SORT_ORDER queryname
    """
}

workflow RerunSortConsensus {
    take:
        tuple_input
        extension
    main:
        SortConsensus(tuple_input, extension)
    emit:
        final_out = SortConsensus.out
}

// Finally, merge the consensus_mapped.bam with the consensus_unmapped.bam to
// retain the UMI group information.
process MergeBam2 {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(bam_aligned), path(bam_unmapped)
        val extension

    output:
        tuple val(sample_id), path("${sample_id}${extension}.ba[m,i]")

    """
    gatk MergeBamAlignment \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_RETAIN RX \
        --ALIGNED_BAM ${bam_aligned} \
        --UNMAPPED_BAM ${bam_unmapped} \
        --OUTPUT ${sample_id}${extension}.bam \
        --REFERENCE_SEQUENCE ${params.ref} \
        --SORT_ORDER coordinate \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --CREATE_INDEX true \
        --CLIP_OVERLAPPING_READS false
    """
}

// Variant calling step using vardict
process VarDict {
    tag "$sample_id"
    
    input:
        tuple val(sample_id), path(bami)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.vcf")

    """
    vardict \
        -G ${params.ref} \
        -f 0.0005 \
        -N ${sample_id} \
        -b ${bami[1]} \
        -c 1 \
        -S 2 \
        -E 3 \
        -g 4 ${params.bed} \
        | $params.teststrandbias \
        | $params.var2vcf > ${sample_id}${extension}.vcf
    """
}

// Annotation step using VEP
process AnnotationVEP {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(vcf)
        val extension
    
    output:
         tuple val(sample_id), file("${sample_id}${extension}.vcf")

    """
    vep \
        -i ${vcf} \
        -o ${sample_id}${extension}.vcf \
        --cache \
        --dir_cache ${params.cache} \
        --offline \
        --force_overwrite \
        --vcf \
        --numbers \
        --refseq \
        --symbol \
        --hgvs \
        --canonical \
        --max_af \
        --fasta ${params.fasta}
    """
}

process BedToIntervalList {

    input:
        path dict
        path bed
        val extension

    output:
        file "${extension}.interval_list"
    
    """
    picard BedToIntervalList \
        --SEQUENCE_DICTIONARY ${dict} \
        --INPUT ${bed} \
        --OUTPUT ${extension}.interval_list
    """
}

// Extraction of read quality metrics before deduplication
process CollectHsMetrics {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), file(bam)
        file bed
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.txt")

    """
    picard CollectHsMetrics \
        --REFERENCE_SEQUENCE ${params.ref} \
        --BAIT_INTERVALS ${bed} \
        --TARGET_INTERVALS ${bed} \
        --INPUT ${bam} \
        --OUTPUT ${sample_id}${extension}.txt
    """
}

workflow RerunCollectHsMetrics {
    take:
        tuple_input
        bed
        extension
    main:
        CollectHsMetrics(tuple_input, bed, extension)
    emit:
        final_out = CollectHsMetrics.out
}

process BCFtools_stats {
    tag "${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(file)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.stats")

    """
    bcftools stats ${file} > ${sample_id}${extension}.stats
    """
}

// Generate a multi-quality control report from collected metrics data (process
// collectmetrics2 output).
process MultiQC {
    tag "${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(file)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}")

    """
    multiqc ${params.outdir}/${sample_id} -o ${sample_id}${extension}
    """
}

workflow {
    Channel.fromFilePairs(params.fastq, checkIfExists:true)
        .filter{ v -> v=~ params.filter_fastq}
        .set{read_pairs_fastq}

    // 1. Preprocess deduplication
    ConvertFastqToSam(read_pairs_fastq, ".1.unmapped")
    ExtractUmis(ConvertFastqToSam.out, ".1.umi_extracted")
    ConvertSamToFastq(ExtractUmis.out, ".1.umi_extracted")
    Fastp(ConvertSamToFastq.out, ".1.umi_extracted.trimmed")
    BWAmem(Fastp.out[0], "-t 10", ".1.umi_extracted.aligned")

    BWAmem.out.join(ExtractUmis.out).set{bams_umis}
    MergeBam(bams_umis, ".1.merged")

    // 2. Process deduplication
    UmiMergeFilt(MergeBam.out, ".2.filtered")
    GroupReads(UmiMergeFilt.out, ".2.umi_grouped")
    CallConsensus(GroupReads.out.nextout, ".2.consensus_unmapped")

    // 3. Post process deduplication
    RerunConvertSamToFastq(CallConsensus.out, ".3.unmapped")
    RerunBWAmem(RerunConvertSamToFastq.out.final_out, "-t 10 -Y", ".3.consensus_mapped")
    SortConsensus(CallConsensus.out, ".3.unmapped")
    RerunSortConsensus(RerunBWAmem.out.final_out, ".3.mapped")

    RerunSortConsensus.out.final_out.join(SortConsensus.out).set{bams_consensus}
    MergeBam2(bams_consensus, ".3.merged")

    // 4. Variant Calling & Annotation
    VarDict(MergeBam2.out, ".4.vardict")
    AnnotationVEP(VarDict.out, ".4.vardict.vep")

    // Quality Controls
    BedToIntervalList(params.dict, params.bed, params.bed.tokenize('.')[0].tokenize('/')[-1])
    CollectHsMetrics(BWAmem.out, BedToIntervalList.out, ".QC.HsMetrics.1")
    RerunCollectHsMetrics(RerunBWAmem.out.final_out, BedToIntervalList.out, ".QC.HsMetrics.3")
    BCFtools_stats(VarDict.out, ".QC.bcftools_stats")
    MultiQC(BCFtools_stats.out, ".QC.multiQC")
}
