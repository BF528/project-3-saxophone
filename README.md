# Project Description

Gene sequencing can be profiled by the traditional microarray or the more advanced and powerful method of RNA seq. RNA-seq can help identify more differentially modulated transcripts of toxicological relevance, splice variants, and non-coding transcripts such as micro RNA (miRNA), long non-coding RNA (lncRNA), pseudogenes. These additional data may be informative for toxicity prediction, mechanistic investigations, or could be the biomarker discovery. 

Wang et al were to study and result in the comparative analysis of gene expression responses profiled by Affymetrix microarray and Illumina RNA-seq in the liver tissues from rats that are exposed to diverse chemicals. 

In this project our goal is the following: 
  1. Align short reads to the rat genome using STAR and quantify expression using a read counting strategy
  2. Perform differential expression analysis of RNA-Seq data with DESeq2
  3. Perform differential expression analysis of pre-normalized microarray expression data using limma
  4. Map between Affymetrix and refSeq identifier systems

# Contributors

Sooyoun Lee leesu@bu.edu - Data Curator 

Jason Rose jjrose@bu.edu - Programmer

Daniel Goldstein djgoldst@bu.edu - Analyst

Sunny Yang yang98@bu.edu - Biologist

# Repository Contents
Data Curator:
1. fastqc.sh

The FastQC module was used to run the quality control checks for all ".fatstq.gz" samples stored in /projectnb/bf528/users/saxophone/project3/ directory. After submitting the fastqc.sh script, the fastq.zip and the fastqc.html for each samples were created. 

2. STAR.sh 

By using the STAR aligner, each of the samples is aligned against the rat genome which can be found in /project/bf528/project_3/reference/rn4_STAR.By changing different "_1.fastq.gz", "_2.fastq.gz"file names, and the outFileNamePrefix names each sample was creating five different files such as Aligned.sortedByCoord.out.bam, Log.final.out, Log .out, Log.progress.out, a,d SJ.out.tab.  

3. multiqc.sh

Summary of combining all the data into a single convenient report. The multiqc.qsub script was submitted and as a result, the multiqc_report.html webpage was created where all the samples are analyzed. 

Programmer:
1. R script for section 3, prepping for DESeq and creating a featurecounts bar plot using the ggplot2 package.

Reads in multiqc feature count output files to prep for DESeq analysis by concatenating all samples counts into one file against matching geneid.

2. R script for section 4, running DESeq and creating histograms and volcano plots using ggplot package.

Reads in the concatenated file from the section 3 R script and sets the gene IDs to row names. Eliminates counts of zero and consturcts the DESeq object. The results are shrunk and volcano plots created for each mode of action. P values are adjusted and counts of significant genes provided at a 0.05 adjp value threshold. Histograms for each mode of action are created.
