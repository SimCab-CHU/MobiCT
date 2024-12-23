
// Generate a multi-quality control report from collected metrics data 
process MultiQC {
    tag "${sample_id}"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(dir)
        val extension

    output:
        file("${sample_id}${extension}")

    """
    multiqc ${dir} -o ${sample_id}${extension} 
    """
}
