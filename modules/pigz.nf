//
process DealWithRef {
    tag "$genome_fasta"
    label 'pigz'
 
    input:
        path(genome_fasta)

    output:
        path "*.fa", emit: genomerefok

    script:
        genomeReady = genome_fasta.baseName.replaceAll(/\.(gz)$/, '') // remove .gz
        genomeReady = genomeReady.replaceAll(/\.(fasta|fa)$/, '') // remove .fasta or .fa
        genomeFa = genomeReady + ".fa"
    """
        

        # DEALING WITH GZIPPED FILES to uncompress if needed
        extension=\$(echo ${genome_fasta} | awk -F . '{print \$NF}')

        if [[ \${extension} == "gz" ]];then
            pigz -dck -p ${task.cpus} ${genome_fasta} > ${genomeFa}
        else
            # copy fa
            ln -s ${genome_fasta} ${genomeFa}
        fi
    """
}