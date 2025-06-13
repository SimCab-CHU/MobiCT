process CreateSequenceDictionary {
    label 'picard'

    input:
        path genome

    output:
        file "*.dict"
    
    """
    picard CreateSequenceDictionary R=${genome} O=${genome}.dict
    """
}

process BedToIntervalList {
    label 'picard'

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
    label 'picard'
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), file(bam)
        file bed
        path genomeref
        path genome_fai
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.txt")

    """
    picard CollectHsMetrics \
        --REFERENCE_SEQUENCE ${genomeref} \
        --BAIT_INTERVALS ${bed} \
        --TARGET_INTERVALS ${bed} \
        --INPUT ${bam} \
        --OUTPUT ${sample_id}${extension}.txt
    """
}
