// Alignement before deduplication: Align quality and adapter trimmed reads to
// the reference genome
process BWAmem {
    tag "$sample_id"

    input:
        tuple val(sample_id), path(fastq)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam")

    script:
        def args = task.ext.args ?: ''

    """
    bwa mem ${args} \
        -t ${task.cpus} \
        -M ${params.ref} \
        ${fastq[0]} \
        ${fastq[1]} \
        | \
        samtools view --threads ${task.cpus-1} -bh -o ${sample_id}${extension}.bam
    """
}
