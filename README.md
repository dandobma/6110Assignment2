# Assignment 2 for BINF6110, Genomic Methods for Bioinformatics

## Matei Dan-Dobre

### 1407506


# Introduction

The ability to quantify genome-wide gene expression has transformed modern biology by enabling systematic investigation of cellular responses to environmental and developmental changes. RNA sequencing (RNA-seq) is now the standard approach for transcriptome profiling due to its high sensitivity, broad dynamic range, and capacity to detect both known and novel transcripts (Wang et al., 2009). In bulk RNA-seq experiments, RNA from multiple cells within a biological condition is sequenced to estimate gene expression levels, allowing statistical comparison between experimental groups. However, RNA-seq data consist of discrete count measurements that exhibit biological and technical variability, necessitating statistical frameworks specifically designed for overdispersed count data. Methods such as DESeq2 model RNA-seq counts using negative binomial distributions and incorporate empirical Bayes shrinkage to provide stable estimates of dispersion and fold change (Love et al., 2014). These approaches enable rigorous identification of differentially expressed genes while controlling for multiple testing.

The budding yeast Saccharomyces cerevisiae serves as an ideal model organism for studying gene regulatory responses to environmental stress and developmental transitions. During biological wine aging, specialized strains of S. cerevisiae form a surface-associated biofilm known as velum. Velum is a membranous structure that develops at the air–liquid interface of fermenting wine and allows yeast cells to persist under oxidative and nutrient-limited conditions. This adaptation is essential for the production of certain aged wines and represents a complex physiological transition involving metabolic reprogramming, stress tolerance, and cell surface remodeling (Alexandre, 2013).

Biofilm formation in yeast is a coordinated developmental process involving adhesion, extracellular matrix production, and global transcriptional shifts. Members of the FLO gene family, particularly FLO11, play central roles in adhesion and biofilm architecture (Verstrepen & Klis, 2006). In addition to adhesion-related genes, successful biofilm development requires extensive metabolic adaptation. As yeast cells transition from fermentative growth to oxidative metabolism, genes involved in respiration, membrane transport, ion homeostasis, and carbohydrate metabolism are frequently differentially regulated (Gasch et al., 2000). Environmental stress responses are also activated under nutrient limitation and oxidative conditions, resulting in coordinated transcriptional programs that enable survival (Gasch et al., 2000). Understanding these regulatory shifts during velum development provides insight into both industrial fermentation processes and broader principles of fungal adaptation.

The dataset analyzed in this study consists of RNA-seq profiles from three developmental stages of velum formation (early, thin, and mature) with biological replicates at each stage. Because gene expression differences across developmental time points may involve hundreds or thousands of genes, systematic statistical comparison is required. Differential expression analysis was performed using DESeq2, which estimates size factors to normalize sequencing depth, models dispersion across biological replicates, and performs hypothesis testing using Wald statistics (Love et al., 2014). Log2 fold-change shrinkage via the apeglm method was applied where appropriate to reduce noise in effect size estimates while preserving biologically meaningful large differences (Zhu et al., 2018). Multiple testing correction was implemented using the Benjamini–Hochberg procedure to control the false discovery rate (Benjamini & Hochberg, 1995), ensuring statistical rigor when evaluating thousands of genes simultaneously.

While differential expression identifies genes whose expression levels change significantly between stages, interpretation of these results requires biological context. Functional annotation assigns genes to known biological categories such as Gene Ontology (GO) terms, protein families, and metabolic pathways. The Gene Ontology provides a structured, hierarchical vocabulary describing gene products in terms of biological processes, molecular functions, and cellular components (Ashburner et al., 2000). Mapping differentially expressed genes to GO categories enables identification of overarching biological themes associated with developmental transitions. Beyond classification, over-representation analysis (ORA) statistically evaluates whether specific GO categories are enriched among differentially expressed genes relative to a genomic background (Subramanian et al., 2005). Enrichment analysis helps distinguish coordinated biological responses from random fluctuations in gene expression.

In the context of velum development, enrichment of processes related to carbohydrate metabolism, membrane transport, stress response, or ion homeostasis would be consistent with known physiological adaptations to oxidative, nutrient-limited environments. Conversely, repression of fermentative pathways may reflect the metabolic shift toward respiratory growth characteristic of surface-associated yeast populations. By integrating differential expression analysis with functional annotation and enrichment testing, it is possible to move from gene-level changes to pathway-level biological interpretation.
The objective of this study is therefore to characterize transcriptional changes across early, thin, and mature stages of velum development using a statistically rigorous RNA-seq workflow. Differential expression analysis was conducted using DESeq2 to identify genes with significant stage-dependent expression changes. Functional annotation was performed to categorize differentially expressed genes by Gene Ontology Biological Process terms, and over-representation analysis was used to identify biological processes significantly enriched at each developmental transition. Together, these approaches provide a comprehensive framework for understanding the molecular mechanisms underlying velum formation and maturation in Saccharomyces cerevisiae.


# Methods

## Computational Environment

All preprocessing and transcript quantification were performed in Windows Subsystem for Linux (WSL) running Ubuntu 22.04 LTS. Downstream statistical and functional analyses were conducted in R (version 4.3.x; R Core Team, 2023). All command-line scripts are provided in the scripts/ directory, and the complete R workflow is documented in scripts/04_deseq2_ora.R.

## Data Acquisition and FASTQ Preparation

Raw RNA-seq data were downloaded from the NCBI Sequence Read Archive (SRA) using accession numbers corresponding to three developmental stages of velum formation (early, thin, and mature; three biological replicates per stage).
Downloads and FASTQ conversion were performed using SRA Toolkit v3.x, specifically the fasterq-dump and prefetch utilities (NCBI SRA Toolkit, 2023). Automated download and compression were handled using:

[`scripts/01_download_sra.sh`](scripts/01_download_sra.sh)

This script attempted direct FASTQ conversion via fasterq-dump and implemented a fallback strategy using prefetch when necessary to ensure robust data retrieval. Resulting FASTQ files were gzip-compressed and stored in data/fastq/.

## Quality Control

Sequencing quality assessment was performed using FastQC v0.11.x (Andrews, 2010), and summary reports were generated using MultiQC v1.x (Ewels et al., 2016). These steps were automated via:

[`scripts/02_qc.sh`](scripts/02_qc.sh)

Quality reports were saved in results/qc/ and inspected to verify overall sequence quality, adapter content, and base quality distribution prior to quantification.

## Reference Preparation

The Saccharomyces cerevisiae reference transcriptome (R64-1-1 cDNA) and corresponding gene annotation (GTF; Ensembl release 111) were downloaded from the Ensembl FTP repository. The transcriptome FASTA file was decompressed and used to construct a Salmon index.

## Salmon Index Construction

Transcriptome indexing was performed using Salmon v1.10.3 (Patro et al., 2017) with a k-mer size of 31:

salmon index \
 -t reference/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa \
 -i reference/salmon_index \
 -k 31
 
Salmon uses a lightweight quasi-mapping approach (pufferfish index) for efficient transcript-level quantification (Patro et al., 2017).

## Transcript Quantification

Transcript abundance was quantified in single-end mode using Salmon v1.10.3 with selective alignment validation and bias correction options enabled (--validateMappings, --gcBias, --seqBias, --posBias). Quantification was automated using:

[`scripts/03_salmon_quant.sh`](scripts/03_salmon_quant.sh)

For each sample, Salmon generated transcript-level abundance estimates in results/salmon/<sample>/quant.sf.

## Transcript-to-Gene Summarization
