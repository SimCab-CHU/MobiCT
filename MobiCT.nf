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

// Check profiles
if (
    workflow.profile.contains('singularity') ||
    workflow.profile.contains('docker') ||
    workflow.profile.contains('conda')
  ) { "executer selected" }
else { exit 1, "No executer selected, please use: -profile docker/singularity/conda (choose one of the 3 possibilities)"}

// Include modules
include {ConvertFastqToSam; MergeBam; MergeBam2; SortConsensus; SortConsensus as RerunSortConsensus; ConvertSamToFastq; ConvertSamToFastq as RerunConvertSamToFastq } from "$baseDir/modules/gatk.nf"
include {CreateSequenceDictionary; BedToIntervalList; CollectHsMetrics; CollectHsMetrics as RerunCollectHsMetrics;} from "$baseDir/modules/picard.nf"
include {BwaIndex; BWAmem; BWAmem as RerunBWAmem;} from "$baseDir/modules/bwamem.nf"
include {BCFtools_stats} from "$baseDir/modules/bcftools.nf"
include {Fastp} from "$baseDir/modules/fastp.nf"
include {ExtractUmis; GroupReads; CallConsensus; } from "$baseDir/modules/fgbio.nf"
include {MultiQC; MultiQC as MultiQC_ALL; } from "$baseDir/modules/multiqc.nf"
include {Faidx; UmiMergeFilt; Sam2bam; Sam2bam as ReSam2bam} from "$baseDir/modules/samtools.nf"
include {DealWithRef} from "$baseDir/modules/pigz.nf"
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


    // htslib-based tools (e.g., samtools, bcftools) Require BGZF-compressed FASTA and  Use .fai index (samtools faidx)
    DealWithRef(genomeref.collect())
    DealWithRef.out.genomerefok.set{genomerefok}

    Faidx(genomerefok).set{genome_fai}

    // Create the genome index
    BwaIndex(genomerefok)

    // Picard/GATK can use gzipped compressed fasta but need a .dict file 
    // Create the genome dictionary
    CreateSequenceDictionary(genomerefok)
        .set{genome_dict}

    // 1. Preprocess deduplication
    ConvertFastqToSam(read_pairs_fastq, ".1.unmapped")
    ExtractUmis(ConvertFastqToSam.out, params.struct_r1, params.struct_r2, ".1.umi_extracted")
    ConvertSamToFastq(ExtractUmis.out, ".1.umi_extracted")
    Fastp(ConvertSamToFastq.out, ".1.umi_extracted.trimmed")
    BWAmem(Fastp.out[0], genomerefok, BwaIndex.out, "", ".1.umi_extracted.aligned")
    Sam2bam(BWAmem.out.tuple_sample_sam)

    Sam2bam.out.join(ExtractUmis.out).set{bams_umis}
    MergeBam(bams_umis, genomerefok, genome_dict, ".1.merged")
    UmiMergeFilt(MergeBam.out, ".1.filtered")

    // 2. Process deduplication
    GroupReads(UmiMergeFilt.out, ".2.umi_grouped")
    CallConsensus(GroupReads.out.nextout, genome_fai, ".2.consensus_unmapped")

    // 3. Post process deduplication
    RerunConvertSamToFastq(CallConsensus.out, ".3.unmapped")
    RerunBWAmem(RerunConvertSamToFastq.out, genomerefok, BwaIndex.out, " -Y", ".3.consensus_mapped")
    SortConsensus(CallConsensus.out, ".3.unmapped")
    ReSam2bam(RerunBWAmem.out.tuple_sample_sam)
    RerunSortConsensus(ReSam2bam.out, ".3.mapped")

    RerunSortConsensus.out.join(SortConsensus.out).set{bams_consensus}
    MergeBam2(bams_consensus, genomerefok, genome_dict, ".3.merged")

    // 4. Variant Calling & Annotation
    VarDict(MergeBam2.out, genomerefok, genome_fai, params.bed,".4.vardict")
    AnnotationVEP(VarDict.out, genomerefok, ".4.vardict.vep")

    // Quality Controls
    BedToIntervalList(genome_dict, params.bed, params.bed.tokenize('.')[0].tokenize('/')[-1])
    CollectHsMetrics(Sam2bam.out, BedToIntervalList.out, genomerefok, genome_fai, ".QC.HsMetrics.1")
    RerunCollectHsMetrics(ReSam2bam.out, BedToIntervalList.out, genomerefok, genome_fai, ".QC.HsMetrics.3")
    BCFtools_stats(VarDict.out, ".QC.bcftools_stats")
    MultiQC(BCFtools_stats.out, ".QC.multiQC")
    MultiQC_ALL(BCFtools_stats.out.collect(), "all.QC.multiQC")
}
