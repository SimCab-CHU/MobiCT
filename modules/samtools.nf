// 
process UmiMergeFilt {
    tag "$sample_id"
    label 'samtools'

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