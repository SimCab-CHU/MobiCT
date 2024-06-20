# ctSOM
Analysis pipeline for somatic diagnosis on ctDNA
# Introduction
ctSOM is a variant analysis pipeline designed for circulating tumor DNA (ctDNA) obtained through non-invasive liquid biopsy, serving diagnostic, prognostic, and therapeutic purposes.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses conda containers making installation trivial and results highly reproducible.

The strength of this pipeline is attributed to the strategic integration of Unique Molecular Identifiers (UMIs). These molecular barcodes effectively mitigate the challenges posed by PCR duplicates and sequencing errors, contributing to the pipeline's enhanced accuracy and sensitivity in identifying genetic variants

<picture>
 <img alt="the workflow" src="https://github.com/Oussamadl/ctSOM/blob/main/pipeline.png">
</picture>


# Pipeline summary
the pipeline can currently perform the following
- Extraction of the UMI from insert reads and generate consensus reads from them (fgbio )
- Adapter and quality trimming (fastp).
- Map Reads to Reference (bwa mem).
- Process BAM file (gatk SortSam, gatk MergeBamAlignment)
- convertion of file format ( gatk SamToFastq, gatk FastqToSam, vcf2tsv)
- quality control (samtools depth, picard collectHsMetrcis, multiQc)
- variant calling (vardict)
- Annotation (VEP)
# add image of pipeline

# Quick Start
1. Install nextflow.
2. creat this conda environment gatk4: conda create -n gatk4 -c conda-forge -c bioconda gatk4.
3. activate the conda environment you create, then istall the following tools:
  - fgbio v2.1.0  : fulcrumgenomics.github.io/fgbio
  - bwa v0.7.17   :bio-bwa.sourceforge.net
  - fastp v0.23.2 :github.com/OpenGene/fastp
  - samtools v1.17 :htslib.org
  - picard
  - vardict
  - VEP
  - vcf2tsvpy
5. Download the reference genome, and extract all files related to it through indexing: hg38.fa.amb, hg38.fa.ann, hg38.fa.bwt, hg38.fa.fai, hg38.fa.pac, hg38.fa.sa
6. Download the dataset needed for VEP use (Home_sapiens.....fa.gz , and the cache: homo_sapiens_vep_......tar.gz)
7. On nextflow.config file, define all input and output paths of your pipeline 
8. Download the pipeline and Run it on your Dataset:
   * Nextflow -log /output_directory/my.log run ctSOM.nf -with-report /output_directory/report.html -with-timeline /output_directory/timeline.html -with-dag /output_directory/flowchart.dot

9. the tsv file output for each sample serve as an input for "annotExtraction.py" python script, which will sort out an xls file containing: gene name, HGVSc (e.g: c.2450T>G), HGVSp (e.g: p.Leu817Arg) and variant allele frequency (VAF).

# Pipeline output
The outcomes of your execution are stored within the directory you specified using the --outdir parameter. The log, HTML, and DOT files, which are provided as optional outputs, offer you valuable insights into the progress of your run:
- my.log : Logs are important for debugging, tracking the execution progress, and identifying any errors that might occur during the execution of the pipeline
- report.html typically includes information about the execution of each process in the workflow, such as input and output files, execution times, and any errors encountered.
- timeline.html provides a visual representation of the execution timeline of processes in the workflow, showing when each process started and finished.
- flowchart.dot: The **.dot** format is a standard format for describing graphs. it visualizes the workflow's structure and dependencies between processes.
