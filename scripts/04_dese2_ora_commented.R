################################################################################
# scripts/04_deseq2_ora.R
#
# Purpose:
#   End-to-end RNA-seq differential expression workflow for yeast:
#     1) Read sample metadata (early / thin / mature)
#     2) Import Salmon quantification results with tximport
#     3) Build DESeq2 dataset and run DE analysis (~ stage)
#     4) Produce key plots:
#          - PCA
#          - Volcano plots (3 contrasts)
#          - Heatmap of top variable genes
#     5) Perform ORA (GO Biological Process) using clusterProfiler
#     6) Save outputs (CSVs, figures, sessionInfo)
#
# Assumptions / Inputs:
#   - You have already run Salmon quantification and have:
#       results/salmon/<sample>/quant.sf
#   - You have created:
#       data/samples.csv
#   - You have downloaded the Ensembl GTF:
#       reference/Saccharomyces_cerevisiae.R64-1-1.111.gtf.gz
#
# Outputs:
#   - results/deseq2/deseq2_<contrast>.csv
#   - figures/fig1_pca.png
#   - figures/fig2_volcano_<contrast>.png
#   - figures/fig3_heatmap_top50var.png
#   - results/ora/ora_go_bp_<contrast>.csv
#   - figures/fig4_ora_go_bp_<contrast>.png
#   - results/sessionInfo_R.txt
################################################################################


################################################################################
# 0) PACKAGE INSTALLATION (RUN ONCE, THEN COMMENT OUT)
################################################################################

# If BiocManager is not installed, install it (needed for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages used in this pipeline
# NOTE: You usually run this once, not every time you run the script.
BiocManager::install(c(
  "tximport",         # import transcript-level quant into R; summarize to gene
  "DESeq2",           # differential expression analysis for count data
  "apeglm",           # shrinkage estimator for log2 fold-changes (coef-based)
  "clusterProfiler",  # enrichment analysis (ORA/GSEA)
  "org.Sc.sgd.db",    # yeast (S. cerevisiae) annotation database
  "AnnotationDbi",    # mapping utilities for OrgDb objects
  "enrichplot"        # plotting functions for enrichment results
))

# Install CRAN packages used in this pipeline
install.packages(c(
  "tidyverse",        # data wrangling + ggplot2
  "pheatmap",         # heatmaps
  "EnhancedVolcano",  # volcano plot helper
  "matrixStats"       # rowVars() used for top-variable gene selection
))


################################################################################
# 1) LOAD LIBRARIES (SUPPRESS STARTUP MESSAGES)
################################################################################

# Load all required packages while suppressing the usual startup text
suppressPackageStartupMessages({
  
  # tidyverse:
  # - dplyr for data manipulation
  # - ggplot2 for plotting
  # - readr for reading/writing CSVs
  library(tidyverse)
  
  # tximport reads Salmon quant.sf and can summarize transcript -> gene
  library(tximport)
  
  # DESeq2 performs differential expression on gene-level count matrix
  library(DESeq2)
  
  # matrixStats provides rowVars(), used to pick top-variable genes for heatmap
  library(matrixStats)
  
  # heatmap plotting
  library(pheatmap)
  
  # volcano plot plotting
  library(EnhancedVolcano)
  
  # enrichment analysis tools
  library(clusterProfiler)
  
  # yeast annotation database
  library(org.Sc.sgd.db)
  
  # annotation mapping helper functions
  library(AnnotationDbi)
  
  # enrichment plotting (dotplot, etc.)
  library(enrichplot)
})


################################################################################
# 2) INPUTS: READ SAMPLE METADATA + DEFINE GROUPS
################################################################################

# Read sample sheet (must exist):
#   - column 'sample' should match your Salmon folder names
#   - column 'stage' must be early/thin/mature
samples <- readr::read_csv("data/samples.csv", show_col_types = FALSE) %>%
  # Ensure stage is treated as an ordered factor (important for reference level)
  # Reference level = first level = "early"
  mutate(stage = factor(stage, levels = c("early", "thin", "mature")))


################################################################################
# 3) INPUTS: LOCATE SALMON QUANT FILES (quant.sf)
################################################################################

# Build a vector of file paths to each sample’s quant.sf file
# The structure assumed:
#   results/salmon/<sample_name>/quant.sf
files <- file.path("results/salmon", samples$sample, "quant.sf")

# Name the vector elements by sample ID so tximport can track sample order
names(files) <- samples$sample

# Stop immediately if any quant.sf files are missing
stopifnot(all(file.exists(files)))


################################################################################
# 4) BUILD tx2gene MAPPING FROM ENSEMBL GTF
#    Why:
#      Salmon outputs transcript-level abundance; DESeq2 expects gene-level counts.
#      tximport needs a transcript->gene mapping (tx2gene) to aggregate properly.
################################################################################

# Path to Ensembl gene annotation (GTF)
# IMPORTANT: The gene_id and transcript_id formats must match Salmon’s transcript IDs
gtf <- "reference/Saccharomyces_cerevisiae.R64-1-1.111.gtf.gz"

# Read the GTF as a tab-delimited file
# - comment = "#" skips header lines
# - col_names = FALSE because GTF doesn't have a header row
gtf_tbl <- readr::read_tsv(
  gtf,
  comment = "#",
  col_names = FALSE,
  show_col_types = FALSE
)

# Assign meaningful column names for standard GTF fields
colnames(gtf_tbl) <- c(
  "seqname",   # chromosome / contig name
  "source",    # source of annotation
  "feature",   # feature type (gene/transcript/exon/etc.)
  "start",     # start coordinate
  "end",       # end coordinate
  "score",     # optional score field
  "strand",    # + or -
  "frame",     # reading frame (for CDS)
  "attribute"  # key-value attributes (gene_id, transcript_id, etc.)
)

# Create transcript -> gene mapping by:
#  - keeping only rows where feature == "transcript"
#  - extracting transcript_id and gene_id from the attribute string
#  - removing rows missing either ID
#  - dropping duplicates
tx2gene <- gtf_tbl %>%
  dplyr::filter(feature == "transcript") %>%   # keep transcript entries only
  dplyr::transmute(
    # Extract transcript_id using a regex group
    transcript_id = stringr::str_match(attribute, 'transcript_id "([^"]+)"')[, 2],
    # Extract gene_id using a regex group
    gene_id       = stringr::str_match(attribute, 'gene_id "([^"]+)"')[, 2]
  ) %>%
  dplyr::filter(!is.na(transcript_id), !is.na(gene_id)) %>%  # drop missing IDs
  dplyr::distinct()                                          # ensure uniqueness


################################################################################
# 5) tximport: IMPORT SALMON + SUMMARIZE TRANSCRIPT -> GENE
################################################################################

# tximport reads each quant.sf and produces:
#   - counts (estimated counts)
#   - abundance (TPM)
#   - length (effective transcript lengths)
#
# Key arguments:
#   type="salmon"  -> tells tximport the file format
#   tx2gene        -> transcript->gene mapping for aggregation
#   ignoreTxVersion=TRUE -> ignores transcript version suffixes if present
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)


################################################################################
# 6) DESeq2 SETUP: CREATE colData + DESeqDataSet
################################################################################

# Build colData (sample metadata) with rownames = sample IDs
# DESeq2 requires rownames in colData to match column names in count matrix
coldata <- samples %>%
  dplyr::select(sample, stage) %>%               # keep only sample + stage
  tibble::column_to_rownames("sample")           # set sample as rownames

# Create DESeqDataSet from tximport output
# design = ~ stage means we model gene expression as a function of stage
dds <- DESeqDataSetFromTximport(
  txi,
  colData = coldata,
  design = ~ stage
)

# Filter low-count genes:
# rowSums(counts(dds)) >= 10 keeps genes with at least 10 total counts across all samples
# This reduces noise and helps DESeq2 estimation
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run the full DESeq2 pipeline:
# - estimate size factors
# - estimate dispersions
# - fit the negative binomial model
# - compute Wald tests for coefficients
dds <- DESeq(dds)

# Print out available coefficient names:
# This helps explain why we have coefficients vs early only, not mature vs thin directly
resultsNames(dds)


################################################################################
# 7) OUTPUT FOLDERS: CREATE DIRECTORIES FOR RESULTS + FIGURES
################################################################################

# Create result directories if they don't already exist
dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)  # DE tables
dir.create("results/ora", recursive = TRUE, showWarnings = FALSE)     # enrichment
dir.create("figures", recursive = TRUE, showWarnings = FALSE)         # figures


################################################################################
# 8) FIGURE 1: PCA ON VST-TRANSFORMED COUNTS
#    Why:
#      PCA shows whether samples cluster by stage and flags outliers/batch effects.
################################################################################

# Variance-stabilizing transformation improves PCA interpretability
# blind=FALSE respects the experimental design when estimating transformation
vsd <- vst(dds, blind = FALSE)

# plotPCA returns PCA coordinates if returnData=TRUE
pca_df <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)

# Percent variance explained by the PCs (stored as attribute)
percentVar <- round(100 * attr(pca_df, "percentVar"))

# Build PCA plot with ggplot
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = stage)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

# Save PCA plot as PNG
ggsave("figures/fig1_pca.png", p_pca, width = 6, height = 5, dpi = 300)


################################################################################
# 9) DIFFERENTIAL EXPRESSION: 3 CONTRASTS + VOLCANO PLOTS
#    Contrasts:
#      1) thin vs early
#      2) mature vs early
#      3) mature vs thin
#
#    Important note on shrinkage:
#      - apeglm shrinkage requires a *coef* that exists in resultsNames(dds)
#      - Because "early" is the reference, DESeq2 provides:
#           stage_thin_vs_early, stage_mature_vs_early
#        but NOT stage_mature_vs_thin (must be computed via contrast directly)
################################################################################

# Define the contrasts explicitly for clarity
contrasts <- list(
  thin_vs_early   = c("stage", "thin", "early"),
  mature_vs_early = c("stage", "mature", "early"),
  mature_vs_thin  = c("stage", "mature", "thin")
)

# Store DE results tables for later ORA step
res_list <- list()

# Map only the contrasts that exist as explicit model coefficients
# (because apeglm can shrink only when coef is available)
coef_map <- c(
  thin_vs_early   = "stage_thin_vs_early",
  mature_vs_early = "stage_mature_vs_early"
)

# Loop through contrasts and compute results + volcano plots
for (nm in names(contrasts)) {
  
  # Get the contrast triple
  con <- contrasts[[nm]]
  
  # Decide whether we can use apeglm shrinkage (coef-based) or not
  if (nm %in% names(coef_map)) {
    
    # CASE A: Shrinkable comparison (vs reference "early")
    # lfcShrink with type="apeglm" reduces noise in log2FC while preserving large effects
    res_shrunk <- lfcShrink(dds, coef = coef_map[[nm]], type = "apeglm")
    
    # Convert results object to a tidy tibble and sort by adjusted p-value
    res_tbl <- as.data.frame(res_shrunk) %>%
      rownames_to_column("gene_id") %>%  # gene IDs become a column
      as_tibble() %>%                    # convert to tibble for tidyverse
      arrange(padj)                      # rank by padj (BH-adjusted p-value)
    
  } else {
    
    # CASE B: mature_vs_thin (no apeglm coef exists)
    # We compute directly using DESeq2 results(..., contrast=...)
    # Note: this is unshrunken log2FC unless you use type="ashr" (optional).
    res <- results(dds, contrast = con)
    
    # Convert to tidy tibble and sort by adjusted p-value
    res_tbl <- as.data.frame(res) %>%
      rownames_to_column("gene_id") %>%
      as_tibble() %>%
      arrange(padj)
  }
  
  # Save DE table for this contrast
  write_csv(res_tbl, file.path("results/deseq2", paste0("deseq2_", nm, ".csv")))
  
  # Store results for later ORA step
  res_list[[nm]] <- res_tbl
  
  # ---- Figure 2: Volcano plot ----
  # We open a PNG device, print the volcano plot, then close device
  png(file.path("figures", paste0("fig2_volcano_", nm, ".png")),
      width = 1800, height = 1400, res = 250)
  
  print(EnhancedVolcano(
    res_tbl,                        # data frame with DE results
    lab = res_tbl$gene_id,          # label points by gene_id (can be crowded)
    x = "log2FoldChange",           # x-axis = log2 fold change
    y = "padj",                     # y-axis = adjusted p-value
    pCutoff = 0.05,                 # significance threshold on padj
    FCcutoff = 1.0,                 # threshold on abs(log2FC) (>= 1)
    title = nm                      # plot title = contrast name
  ))
  
  dev.off()  # close PNG device
}


################################################################################
# 10) FIGURE 3: HEATMAP OF TOP 50 MOST VARIABLE GENES (VST)
#     Why:
#       Highlights global expression patterns across samples and stages.
################################################################################

# Identify the 50 genes with highest variance across samples using VST expression
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

# Extract expression matrix for those genes
mat <- assay(vsd)[topVarGenes, , drop = FALSE]

# Mean-center each row (gene) so heatmap shows relative expression differences
mat <- mat - rowMeans(mat)

# Save heatmap as PNG
png("figures/fig3_heatmap_top50var.png", width = 1800, height = 1600, res = 250)

# Plot heatmap with column annotation for stage
pheatmap(mat, annotation_col = as.data.frame(coldata["stage"]))

dev.off()  # close PNG device

################################################################################
# 10.5) Functional annotation (gene-level)
#       This is NOT enrichment.
#       It maps gene IDs -> gene names/descriptions + GO/PFAM/pathway annotations.
################################################################################

# Create output folder for annotated DE tables
dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)

# Choose the keytype that matches your gene_id column
# Most yeast gene IDs look like YAL001C etc. (ORF), so ORF is usually correct.
keytype_use <- "ORF"

# (Optional debug check — comment out once confirmed)
# head(res_list[[1]]$gene_id, 20)

# Function: annotate a vector of gene IDs
annotate_genes <- function(gene_ids, keytype_use = "ORF") {
  
  # Ensure IDs are unique (prevents unnecessary repeated lookups)
  gene_ids <- unique(gene_ids)
  
  # OPTIONAL: remove "-A", "-B" suffixes for annotation lookup only
  # (helps mapping of alternative ORFs like YPR036W-A)
  # gene_ids_lookup <- sub("-[A-Z]$", "", gene_ids)
  
  # ---- 1) Basic 1:1-ish annotation (safe/small) ----
  basic_cols <- c("GENENAME", "COMMON", "DESCRIPTION", "SGD", "UNIPROT")
  
  anno_basic <- AnnotationDbi::select(
    x       = org.Sc.sgd.db,
    keys    = gene_ids,
    keytype = keytype_use,
    columns = basic_cols
  ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(gene_id = !!keytype_use) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      GENENAME    = dplyr::first(na.omit(GENENAME)),
      COMMON      = dplyr::first(na.omit(COMMON)),
      DESCRIPTION = dplyr::first(na.omit(DESCRIPTION)),
      SGD         = dplyr::first(na.omit(SGD)),
      UNIPROT     = dplyr::first(na.omit(UNIPROT)),
      .groups = "drop"
    )
  
  # ---- 2) GO + ONTOLOGY (many-to-one; collapsed lists) ----
  anno_go <- AnnotationDbi::select(
    x       = org.Sc.sgd.db,
    keys    = gene_ids,
    keytype = keytype_use,
    columns = c("GO", "ONTOLOGY")
  ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(gene_id = !!keytype_use) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      GO       = paste(sort(unique(na.omit(GO))), collapse = ";"),
      ONTOLOGY = paste(sort(unique(na.omit(ONTOLOGY))), collapse = ";"),
      .groups = "drop"
    )
  
  # ---- 3) PATH, PFAM, INTERPRO (many-to-one; collapsed lists) ----
  anno_other <- AnnotationDbi::select(
    x       = org.Sc.sgd.db,
    keys    = gene_ids,
    keytype = keytype_use,
    columns = c("PATH", "PFAM", "INTERPRO")
  ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(gene_id = !!keytype_use) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      PATH     = paste(sort(unique(na.omit(PATH))), collapse = ";"),
      PFAM     = paste(sort(unique(na.omit(PFAM))), collapse = ";"),
      INTERPRO = paste(sort(unique(na.omit(INTERPRO))), collapse = ";"),
      .groups = "drop"
    )
  
  # ---- 4) Merge into one annotation table ----
  anno_basic %>%
    dplyr::left_join(anno_go, by = "gene_id") %>%
    dplyr::left_join(anno_other, by = "gene_id")
}

# Annotate + save an annotated DE table for each contrast
for (nm in names(res_list)) {
  
  res_tbl <- res_list[[nm]]
  
  sig_genes <- res_tbl %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  # Skip if there are no significant genes under your thresholds
  if (length(sig_genes) == 0) next
  
  ann_sig <- annotate_genes(sig_genes, keytype_use = keytype_use)
  
  res_annotated <- res_tbl %>%
    dplyr::left_join(ann_sig, by = "gene_id") %>%
    dplyr::arrange(padj)
  
  readr::write_csv(
    res_annotated,
    file.path("results/annotation", paste0("deseq2_", nm, "_annotated.csv"))
  )
}

head(res_list[[1]]$gene_id, 20)
################################################################################
# 11) ORA: GO BIOLOGICAL PROCESS ENRICHMENT (clusterProfiler)
#     Why:
#       Tests whether significant DE genes are overrepresented in GO BP categories.
#
#     Key choices:
#       - Input gene list = significant genes (padj < 0.05 and |log2FC| >= 1)
#       - Multiple testing correction = BH
#       - ont="BP" focuses on Biological Process
#
#     Note about gene IDs:
#       - org.Sc.sgd.db may or may not support ENSEMBL keys depending on version.
#       - We detect whether ENSEMBL keyType is available; otherwise fall back to ORF.
#       - We set readable=FALSE to avoid symbol conversion errors in some installs.
################################################################################

# List supported key types for this OrgDb
kt <- keytypes(org.Sc.sgd.db)

# Choose key type to match your DE gene_id format
# If org.Sc.sgd.db supports ENSEMBL, use that; otherwise use ORF.
keytype_use <- if ("ENSEMBL" %in% kt) "ENSEMBL" else "ORF"

# Loop through each contrast result table stored in res_list
for (nm in names(res_list)) {
  
  # Pull the DE table for this contrast
  res_tbl <- res_list[[nm]]
  
  # Define the significant gene set:
  # - padj < 0.05 (BH adjusted)
  # - abs(log2FoldChange) >= 1 (at least 2-fold change)
  sig <- res_tbl %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    pull(gene_id) %>%
    unique()
  
  # Skip ORA if too few genes (clusterProfiler becomes unstable / uninformative)
  if (length(sig) < 10) next
  
  # Run GO enrichment (ORA) for Biological Process terms
  ego <- enrichGO(
    gene = sig,               # gene list for enrichment
    OrgDb = org.Sc.sgd.db,    # organism database
    keyType = keytype_use,    # how to interpret the gene IDs
    ont = "BP",               # Biological Process GO terms
    pAdjustMethod = "BH",     # multiple testing correction
    readable = FALSE          # do NOT attempt to convert IDs to symbols
  )
  
  # Skip if enrichment returns no terms
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next
  
  # Save ORA table to CSV
  write_csv(as.data.frame(ego), file.path("results/ora", paste0("ora_go_bp_", nm, ".csv")))
  
  # Plot top enriched categories with a dotplot
  p <- dotplot(ego, showCategory = 15) + ggtitle(paste0("GO BP ORA: ", nm))
  
  # Save ORA dotplot figure
  ggsave(file.path("figures", paste0("fig4_ora_go_bp_", nm, ".png")),
         p, width = 8, height = 6, dpi = 300)
}


################################################################################
# 12) REPRODUCIBILITY: SAVE SESSION INFO
#     Why:
#       Captures R version and package versions used to produce results.
################################################################################

# Write sessionInfo() to a file for reproducibility
writeLines(capture.output(sessionInfo()), "results/sessionInfo_R.txt")





keytypes(org.Sc.sgd.db)
