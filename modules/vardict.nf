// Variant calling step using vardict
process VarDict {
    tag "$sample_id"
    label 'verdict'

    input:
        tuple val(sample_id), path(bami)
        path genomeref
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.vcf")

    """
    vardict \
        -G ${genomeref} \
        -f 0.0005 \
        -N ${sample_id} \
        -b ${bami[1]} \
        -c 1 \
        -S 2 \
        -E 3 \
        -g 4 ${params.bed} \
        | $params.teststrandbias \
        | $params.var2vcf > ${sample_id}${extension}.vcf
    """
}
