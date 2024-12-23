
// Variant calling step using vardict
process VarDict {
    tag "$sample_id"
    
    input:
        tuple val(sample_id), path(bami)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.vcf")

    """
    vardict-java \
        -G ${params.ref} \
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
