/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_STATS                    } from '../modules/nf-core/bcftools/stats/main'
include { BWA_MEM                           } from '../modules/nf-core/bwa/mem/main'
include { BWA_MEM as rerunBWA_MEM           } from '../modules/nf-core/bwa/mem/main'
include { ENSEMBLVEP_VEP                    } from '../modules/nf-core/ensemblvep/vep/main'
include { FASTP                             } from '../modules/nf-core/fastp/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from '../modules/nf-core/fgbio/callmolecularconsensusreads/main'
include { FGBIO_GROUPREADSBYUMI             } from '../modules/nf-core/fgbio/groupreadsbyumi/main'
include { GATK4_FASTQTOSAM                  } from '../modules/nf-core/gatk4/fastqtosam/main'
include { GATK4_SAMTOFASTQ                  } from '../modules/nf-core/gatk4/samtofastq/main'
include { GATK4_SAMTOFASTQ as rerunGATK4_SAMTOFASTQ } from '../modules/nf-core/gatk4/samtofastq/main'
include { GATK4_MERGEBAMALIGNMENT           } from '../modules/nf-core/gatk4/mergebamalignment/main'
include { GATK4_MERGEBAMALIGNMENT as rerunGATK4_MERGEBAMALIGNMENT } from '../modules/nf-core/gatk4/mergebamalignment/main'
include { PICARD_BEDTOINTERVALLIST          } from '../modules/nf-core/picard/bedtointervallist/main'
include { PICARD_COLLECTHSMETRICS           } from '../modules/nf-core/picard/collecthsmetrics/main'
include { PICARD_SORTSAM                    } from '../modules/nf-core/picard/sortsam/main'
include { PICARD_SORTSAM as rerunPICARD_SORTSAM } from '../modules/nf-core/picard/sortsam/main'
include { SAMTOOLS_INDEX                    } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW                     } from '../modules/nf-core/samtools/view/main'
include { VARDICTJAVA                       } from '../modules/nf-core/vardictjava/main'
include { FGBIO_EXTRACTUMISFROMBAM          } from '../modules/local/fgbio/extractumisfrombam/main'

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mobict_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MOBICT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    GATK4_FASTQTOSAM(ch_samplesheet)
    ch_versions = ch_versions.mix(GATK4_FASTQTOSAM.out.versions)

    FGBIO_EXTRACTUMISFROMBAM(GATK4_FASTQTOSAM.out.bam, params.read_structure, params.molecular_index_tags)
    ch_versions = ch_versions.mix(FGBIO_EXTRACTUMISFROMBAM.out.versions)

    GATK4_SAMTOFASTQ(FGBIO_EXTRACTUMISFROMBAM.out.bam)
    ch_versions = ch_versions.mix(GATK4_SAMTOFASTQ.out.versions)

    FASTP(GATK4_SAMTOFASTQ.out.fastq, [], false, false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    BWA_MEM(FASTP.out.reads,["idx", params.index_alignment], ["fasta", params.fasta], false)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    BWA_MEM.out.bam.join(FGBIO_EXTRACTUMISFROMBAM.out.bam).set{bams_umis}
    GATK4_MERGEBAMALIGNMENT(bams_umis, ["fasta", params.fasta], ["dict", params.dict])
    ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT.out.versions)
    
    // NEEDS TO BE IMPROVED - bai expected but cannot be generated due to srt order by 'queryname'
    // As bai is not necessecary, use join with another to avoir emty data expected (and mandatory) by SAMTOOLS_VIEW
    GATK4_MERGEBAMALIGNMENT.out.bam.join( BWA_MEM.out.bam ).set{bam_merged}
    ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT.out.versions)
    SAMTOOLS_VIEW(bam_merged, ["fasta", params.fasta], [])
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    FGBIO_GROUPREADSBYUMI(SAMTOOLS_VIEW.out.bam, "adjacency")
    ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_GROUPREADSBYUMI.out.histogram)

    FGBIO_CALLMOLECULARCONSENSUSREADS(FGBIO_GROUPREADSBYUMI.out.bam, 2, 10)
    ch_versions = ch_versions.mix(FGBIO_CALLMOLECULARCONSENSUSREADS.out.versions)

    rerunGATK4_SAMTOFASTQ(FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam)
    ch_versions = ch_versions.mix(rerunGATK4_SAMTOFASTQ.out.versions)
    rerunBWA_MEM(rerunGATK4_SAMTOFASTQ.out.fastq, ["idx", params.index_alignment], ["fasta", params.fasta], false)
    ch_versions = ch_versions.mix(rerunBWA_MEM.out.versions)

    PICARD_SORTSAM(FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam, "queryname")
    ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions)
    rerunPICARD_SORTSAM(rerunBWA_MEM.out.bam, "queryname")
    ch_versions = ch_versions.mix(rerunPICARD_SORTSAM.out.versions)

    rerunPICARD_SORTSAM.out.bam.join(PICARD_SORTSAM.out.bam).set{ bams_sorted }
    rerunGATK4_MERGEBAMALIGNMENT(bams_sorted, ["fasta", params.fasta], ["dict", params.dict])
    ch_versions = ch_versions.mix(rerunGATK4_MERGEBAMALIGNMENT.out.versions)

    SAMTOOLS_INDEX(rerunGATK4_MERGEBAMALIGNMENT.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    rerunGATK4_MERGEBAMALIGNMENT.out.bam.join(SAMTOOLS_INDEX.out.bai).combine([params.bed]).set{ bams_aligned }
    VARDICTJAVA(bams_aligned, ["fasta", params.fasta], ["fai", params.fai])
    ch_versions = ch_versions.mix(VARDICTJAVA.out.versions)

    VARDICTJAVA.out.vcf.join(BWA_MEM.out.bam ).set{vcf}
    ENSEMBLVEP_VEP(vcf, params.assembly, params.species, params.cache_version, params.cache_dir, ["fasta", params.fasta], params.index_alignment)
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
