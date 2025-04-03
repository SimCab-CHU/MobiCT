// Adapter and quality trimming
process Fastp {
    tag "$sample_id"
    label 'fastp'

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.{json,html}'

    input:
        tuple val(sample_id), path(fastq)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.R[1,2].fq")
        file "${sample_id}.QC.fastp.json"
        file "${sample_id}.QC.fastp.html"

    script:
    """
    fastp \
        --thread ${task.cpus} \
        -i ${fastq[0]} \
        -o ${sample_id}${extension}.R1.fq \
        -I ${fastq[1]} \
        -O ${sample_id}${extension}.R2.fq \
        -g -W 5 -q 20 -u 40 -x -3 -l 75 -c \
        -j ${sample_id}.QC.fastp.json \
        -h ${sample_id}.QC.fastp.html \
        -w 12
    """
}