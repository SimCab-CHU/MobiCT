// Variant calling step using vardict
process VarDict {
    tag "$sample_id"
    label 'vardict'

    input:
        tuple val(sample_id), path(bami)
        path genomeref
        path genomeref_fai
        path bed
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.vcf")

    """
    vardict-java \
        -G ${genomeref} \
        -f 0.0005 \
        -N ${sample_id} \
        -b ${bami[1]} \
        -c 1 \
        -S 2 \
        -E 3 \
        -g 4 ${bed} \
        | teststrandbias.R \
        | var2vcf_valid.pl > ${sample_id}${extension}.vcf
    """
}
