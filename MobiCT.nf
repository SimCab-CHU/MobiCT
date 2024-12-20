// PROCESSES

include { ConvertFastqToSam } from "./includes/gatk4.nf"
include { ConvertSamToFastq } from "./includes/gatk4.nf"
include { MergeBam } from "./includes/gatk4.nf"
include { ConvertSamToFastq as RerunConvertSamToFastq } from "./includes/gatk4.nf"
include { SortConsensus } from "./includes/gatk4.nf"
include { SortConsensus as RerunSortConsensus } from "./includes/gatk4.nf"
include { MergeBam2 } from "./includes/gatk4.nf"
include { BedToIntervalList } from "./includes/gatk4.nf"
include { CollectHsMetrics } from "./includes/gatk4.nf"
include { CollectHsMetrics as RerunCollectHsMetrics } from "./includes/gatk4.nf"

include { ExtractUmis } from "./includes/fgbio.nf"
include { GroupReads } from "./includes/fgbio.nf"
include { CallConsensus } from "./includes/fgbio.nf"

include { BWAmem } from "./includes/bwa.nf"
include { BWAmem as RerunBWAmem } from "./includes/bwa.nf"

include { Fastp } from "./includes/fastp.nf"

include { UmiMergeFilt } from "./includes/samtools.nf"

include { VarDict } from "./includes/vardict.nf"

include { AnnotationVEP } from "./includes/vep.nf"

include { BCFtools_stats } from "./includes/bcftools.nf"

include { MultiQC } from "./includes/multiqc.nf"
include { MultiQC as MultiQC_ALL } from "./includes/multiqc.nf"

workflow {
    Channel.fromFilePairs(params.fastq, checkIfExists:true)
        .filter{ v -> v=~ params.filter_fastq}
        .set{read_pairs_fastq}

    // 1. Preprocess deduplication
    ConvertFastqToSam(read_pairs_fastq, ".1.unmapped")
    ExtractUmis(ConvertFastqToSam.out, params.struct_r1, params.struct_r2, ".1.umi_extracted")
    ConvertSamToFastq(ExtractUmis.out, ".1.umi_extracted")
    Fastp(ConvertSamToFastq.out, ".1.umi_extracted.trimmed")
    BWAmem(Fastp.out[0], "-t 10", ".1.umi_extracted.aligned")

    BWAmem.out.join(ExtractUmis.out).set{bams_umis}
    MergeBam(bams_umis, ".1.merged")
    UmiMergeFilt(MergeBam.out, ".1.filtered")

    // 2. Process deduplication
    GroupReads(UmiMergeFilt.out, ".2.umi_grouped")
    CallConsensus(GroupReads.out.nextout, ".2.consensus_unmapped")

    // 3. Post process deduplication
    RerunConvertSamToFastq(CallConsensus.out, ".3.unmapped")
    RerunBWAmem(RerunConvertSamToFastq.out, "-t 10 -Y", ".3.consensus_mapped")
    SortConsensus(CallConsensus.out, ".3.unmapped")
    RerunSortConsensus(RerunBWAmem.out, ".3.mapped")

    RerunSortConsensus.out.join(SortConsensus.out).set{bams_consensus}
    MergeBam2(bams_consensus, ".3.merged")

    // 4. Variant Calling & Annotation
    VarDict(MergeBam2.out, ".4.vardict")
    AnnotationVEP(VarDict.out, ".4.vardict.vep")

    // Quality Controls
    BedToIntervalList(params.dict, params.bed, params.bed.tokenize('.')[0].tokenize('/')[-1])
    CollectHsMetrics(BWAmem.out, BedToIntervalList.out, ".QC.HsMetrics.1")
    RerunCollectHsMetrics(RerunBWAmem.out, BedToIntervalList.out, ".QC.HsMetrics.3")
    BCFtools_stats(VarDict.out, ".QC.bcftools_stats")

    BCFtools_stats.out\
        .map{ it -> [it[0], "$params.outdir/${it[0]}"]}\
        .groupTuple()\
        .set{ multiqc_tuples }

    MultiQC(multiqc_tuples, ".QC.multiQC")

    MultiQC.out.map{ it -> ["all", params.outdir]}.first().set{ all_qc }    
    MultiQC_ALL(all_qc, ".QC.multiQC")
}
