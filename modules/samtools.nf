//
process Faidx {
    tag "$genome_fasta"
    label 'samtools'

    input:
        path(genome_fasta)

    output:
        path("*")

    """
    samtools faidx $genome_fasta
    """
}

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

process Sam2bam {
    label 'samtools'
    tag "$sample"

    input:
        tuple val(sample), path(sam)

    output:
        tuple val(sample), path ("*.bam"), emit: tuple_sample_bam

    script:

        """
            samtools view -@ ${task.cpus} ${sam} -bh -o ${sam.baseName}.bam 
        """

}