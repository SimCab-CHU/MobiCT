
process UmiMergeFilt {
    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam")

    """
    samtools view \
        --threads ${task.cpus -1} \
        -f2 \
        -bh ${bam} \
        -o ${sample_id}${extension}.bam
    """
}
