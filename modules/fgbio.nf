// Extraction of UMIs from the insert reads
process ExtractUmis {
    tag "$sample_id"
    label 'fgbio'

    input:
        tuple val(sample_id), path(bam_file)
        val struct_r1
        val struct_r2
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam")

    """
    fgbio ExtractUmisFromBam \
        -i ${bam_file} \
        -o ${sample_id}${extension}.bam \
        -r ${struct_r1} ${struct_r2} \
        -t RX \
        -a true
    """
}

// Identify and group reads originating from the same source molecule
// The user can control how many errors/mismatches are allowed in the UMI sequence when assigning source molecules (--edits=n).
process GroupReads {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.txt'

    input:
        tuple val(sample_id), path(bam)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam"), emit: nextout
        file "${sample_id}.QC.family_size_counts.txt"

    """
    fgbio GroupReadsByUmi \
        -i ${bam} \
        -o ${sample_id}${extension}.bam \
        --strategy=adjacency \
        --edits=1 \
        -t RX \
        -f ${sample_id}.QC.family_size_counts.txt
    """
}

// Generate consensus reads
// Calculate the consensus sequence, Reads that occur as singletons are
// discarded by default but this can be changed by setting the –min-reads flag
// to 1, in so doing the single read will be considered the consensus.
process CallConsensus {
    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam)
        val extension

    output:
        tuple val(sample_id), file("${sample_id}${extension}.bam")

    """
    fgbio CallMolecularConsensusReads \
        -i ${bam} \
        -o ${sample_id}${extension}.bam \
        --error-rate-post-umi 40 \
        --error-rate-pre-umi 45 \
        --output-per-base-tags false \
        --min-reads 2 \
        --max-reads 50 \
        --read-name-prefix='consensus'
    """
}