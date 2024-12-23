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
    gatk CollectHsMetrics \
        --REFERENCE_SEQUENCE ${params.ref} \
        --BAIT_INTERVALS ${bed} \
        --TARGET_INTERVALS ${bed} \
        --INPUT ${bam} \
        --OUTPUT ${sample_id}${extension}.txt
    """
}
