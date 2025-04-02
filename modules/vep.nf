// Annotation step using VEP
process AnnotationVEP {
    tag "$sample_id"
    label 'vep'

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(vcf)
        val extension
    
    output:
         tuple val(sample_id), file("${sample_id}${extension}.vcf")

    """
    vep \
        -i ${vcf} \
        -o ${sample_id}${extension}.vcf \
        --cache \
        --dir_cache ${params.cache} \
        --offline \
        --force_overwrite \
        --vcf \
        --numbers \
        --refseq \
        --symbol \
        --hgvs \
        --canonical \
        --max_af \
        --fasta ${params.fasta}
    """
}