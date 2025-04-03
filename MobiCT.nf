/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Copyright (C) 2024

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
authors = 'Simon Cabello-Aguilar', 'Charles Van Gothem', 'Jean-Charles Delmas', 'Oussama Bourbia'
copyright = 'Copyright (C) 2024'
license = 'GNU General Public License'
version = '1.0.0'
email = 's-cabelloaguilar@chu-montpellier.fr'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Include modules
include {ConvertFastqToSam; MergeBam; MergeBam2; SortConsensus; SortConsensus as RerunSortConsensus; ConvertSamToFastq; ConvertSamToFastq as RerunConvertSamToFastq } from "$baseDir/modules/gatk.nf"
include {BedToIntervalList; CollectHsMetrics; CollectHsMetrics as RerunCollectHsMetrics;} from "$baseDir/modules/picard.nf"
include {BWAmem; BWAmem as RerunBWAmem;} from "$baseDir/modules/bwamem.nf"
include {BCFtools_stats} from "$baseDir/modules/bcftools.nf"
include {Fastp} from "$baseDir/modules/fastp.nf"
include {ExtractUmis; GroupReads; CallConsensus; } from "$baseDir/modules/fgbio.nf"
include {MultiQC; MultiQC_ALL; } from "$baseDir/modules/multiqc.nf"
include {UmiMergeFilt} from "$baseDir/modules/samtools.nf"
include {VarDict} from "$baseDir/modules/vardict.nf"
include {AnnotationVEP} from "$baseDir/modules/vep.nf"


// Define the workflow
workflow {
    Channel.fromFilePairs(params.fastq, checkIfExists:true)
        .filter{ v -> v=~ params.filter_fastq}
        .set{read_pairs_fastq}
    Channel.fromPath(params.ref, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.ref}!\n" }
            .set{genomeref}

    // 1. Preprocess deduplication
    ConvertFastqToSam(read_pairs_fastq, ".1.unmapped")
    ExtractUmis(ConvertFastqToSam.out, params.struct_r1, params.struct_r2, ".1.umi_extracted")
    ConvertSamToFastq(ExtractUmis.out, ".1.umi_extracted")
    Fastp(ConvertSamToFastq.out, ".1.umi_extracted.trimmed")
    BWAmem(Fastp.out[0], genomeref, "-t 10", ".1.umi_extracted.aligned")

    BWAmem.out.join(ExtractUmis.out).set{bams_umis}
    MergeBam(bams_umis, genomeref, ".1.merged")
    UmiMergeFilt(MergeBam.out, ".1.filtered")

    // 2. Process deduplication
    GroupReads(UmiMergeFilt.out, ".2.umi_grouped")
    CallConsensus(GroupReads.out.nextout, ".2.consensus_unmapped")

    // 3. Post process deduplication
    RerunConvertSamToFastq(CallConsensus.out, ".3.unmapped")
    RerunBWAmem(RerunConvertSamToFastq.out, genomeref, "-t 10 -Y", ".3.consensus_mapped")
    SortConsensus(CallConsensus.out, ".3.unmapped")
    RerunSortConsensus(RerunBWAmem.out, ".3.mapped")

    RerunSortConsensus.out.join(SortConsensus.out).set{bams_consensus}
    MergeBam2(bams_consensus, genomeref, ".3.merged")

    // 4. Variant Calling & Annotation
    VarDict(MergeBam2.out, genomeref, ".4.vardict")
    AnnotationVEP(VarDict.out, ".4.vardict.vep")

    // Quality Controls
    BedToIntervalList(params.dict, params.bed, params.bed.tokenize('.')[0].tokenize('/')[-1])
    CollectHsMetrics(BWAmem.out, BedToIntervalList.out, genomeref, ".QC.HsMetrics.1")
    RerunCollectHsMetrics(RerunBWAmem.out, BedToIntervalList.out, genomeref, ".QC.HsMetrics.3")
    BCFtools_stats(VarDict.out, ".QC.bcftools_stats")
    MultiQC(BCFtools_stats.out, ".QC.multiQC")

    Channel.empty()
        .mix( MultiQC.out )
        .map { sample, files -> files }
        .collect()
        .set { log_files }
    MultiQC_ALL(log_files, "all.QC.multiQC")
}
