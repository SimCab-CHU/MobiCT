// Annotation step using VEP
process AnnotationVEP {
    tag "$sample_id"
    label 'vep'

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(vcf)
        path genomeref
        val extension
    
    output:
         tuple val(sample_id), file("${sample_id}${extension}.vcf")

    script:
        def config = params.cache?.trim() ? "--offline --cache --dir_cache ${params.cache} --max_af " : "--database" 
    """
    vep \
        ${config} \
        --fork ${task.cpus} \
        -i ${vcf} \
        -o ${sample_id}${extension}.vcf \
        --force_overwrite \
        --vcf \
        --numbers \
        --refseq \
        --symbol \
        --hgvs \
        --canonical \
        --fasta ${genomeref}
    """
}