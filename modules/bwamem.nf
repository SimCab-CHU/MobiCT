process BwaIndex {
    label 'bwamem'
    tag "$genome_fasta"

    input:
        path(genome_fasta)

    output:
        path("*")

    script:
        """
        bwa index $genome_fasta -p ${genome_fasta.baseName}
        """
}

// Alignement before deduplication: Align quality and adapter trimmed reads to
// the reference genome
process BWAmem {
    tag "$sample_id"
    label 'bwamem'

    input:
        tuple val(sample_id), path(fastq)
        path genomeref
        path bwa_index_files
        val opt_bwa
        val extension
    
    output:
        tuple val(sample_id), file("${sample_id}${extension}.sam"), emit: tuple_sample_sam
        path "*.log",  emit: bwamem_summary

    """
    bwa mem \
        -t ${task.cpus} \
        ${opt_bwa} \
        -M \
        ${genomeref.baseName} \
        ${fastq[0]} \
        ${fastq[1]} \
        > ${sample_id}${extension}.sam 2> ${sample_id}${extension}.log 
    """
}