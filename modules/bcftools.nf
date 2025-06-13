process BCFtools_stats {
    tag "${sample_id}"
    label 'bcftools'

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(file)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.stats")

    """
    bcftools stats --threads ${task.cpus} ${file} > ${sample_id}${extension}.stats
    """
}
