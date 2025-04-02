// Generate a multi-quality control report from collected metrics data (process
// collectmetrics2 output).
process MultiQC {
    tag "${sample_id}"
    label 'multiqc'

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(file)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}")

    """
    multiqc ${params.outdir}/${sample_id} -o ${sample_id}${extension}
    """
}

process MultiQC_ALL {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(file)
        val extension

    output:
        file("${extension}")

    """
    multiqc ${params.outdir} -o ${extension}
    """
}