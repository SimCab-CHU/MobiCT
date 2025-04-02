// Alignement before deduplication: Align quality and adapter trimmed reads to
// the reference genome
process BWAmem {
    tag "$sample_id"
    label 'bwamem'

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