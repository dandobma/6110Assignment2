####scripts/04_deseq2_ora.R####
# Matei Dan-Dobre
# 1407506


####0) PACKAGE INSTALLATION####
# (RUN ONCE, THEN COMMENT OUT)

#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}
#
#BiocManager::install(c(
#  "tximport",
#  "DESeq2",
#  "apeglm",
#  "clusterProfiler",
#  "org.Sc.sgd.db",
#  "AnnotationDbi",
#  "enrichplot"
#))
#
#install.packages(c(
#  "tidyverse",
#  "pheatmap",
#  "EnhancedVolcano",
#  "matrixStats"
#))



####1) LOAD LIBRARIES (SUPPRESS STARTUP MESSAGES)####

suppressPackageStartupMessages({
  
  library(tidyverse)        # data wrangling + ggplot2
  library(tximport)         # import Salmon quantifications
  library(DESeq2)           # differential expression analysis
  library(matrixStats)      # rowVars() for top-variable genes
  library(pheatmap)         # heatmap plotting
  library(EnhancedVolcano)  # volcano plot helper
  
  library(clusterProfiler)  # ORA / enrichment analysis
  library(org.Sc.sgd.db)    # S. cerevisiae annotation database
  library(AnnotationDbi)    # select() mapping helper used with org.Sc.sgd.db
  library(enrichplot)       # plotting enrichment results (dotplot)
})


####2) INPUTS: READ SAMPLE METADATA + DEFINE GROUPS####
# Read sample sheet and ensure stage factor levels are ordered so "early" is the reference.

samples <- readr::read_csv("data/samples.csv", show_col_types = FALSE) %>%
  mutate(stage = factor(stage, levels = c("early", "thin", "mature")))


####3) INPUTS: LOCATE SALMON QUANT FILES (quant.sf)#####
# Build expected paths to Salmon quant.sf outputs and verify they all exist.

files <- file.path("results/salmon", samples$sample, "quant.sf")  # one path per sample
names(files) <- samples$sample                                    # name vector for tximport sample tracking
stopifnot(all(file.exists(files)))                                # stop early if any file is missing


####4) BUILD tx2gene MAPPING FROM ENSEMBL GTF####
# Parse the Ensembl GTF to construct transcript_id -> gene_id mapping for tximport.

gtf <- "reference/Saccharomyces_cerevisiae.R64-1-1.111.gtf.gz"

gtf_tbl <- readr::read_tsv(
  gtf,
  comment = "#",          # ignore GTF header/comment lines
  col_names = FALSE,      # GTF has no header row
  show_col_types = FALSE
)

# Assign standard GTF column names
colnames(gtf_tbl) <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

# Keep transcript features only, then extract transcript_id and gene_id from the attribute field.
tx2gene <- gtf_tbl %>%
  dplyr::filter(feature == "transcript") %>%
  dplyr::transmute(
    transcript_id = stringr::str_match(attribute, 'transcript_id "([^"]+)"')[, 2],
    gene_id       = stringr::str_match(attribute, 'gene_id "([^"]+)"')[, 2]
  ) %>%
  dplyr::filter(!is.na(transcript_id), !is.na(gene_id)) %>%  # discard incomplete entries
  dplyr::distinct()                                          # remove duplicates


####5) tximport: IMPORT SALMON + SUMMARIZE TRANSCRIPT -> GENE####
# Import Salmon quantifications and aggregate transcript-level counts to gene-level counts.

txi <- tximport(
  files,
  type = "salmon",        # tells tximport what format the input is
  tx2gene = tx2gene,      # transcript -> gene mapping for summarization
  ignoreTxVersion = TRUE  # ignore transcript version suffixes if present
)


####6) DESeq2 SETUP: CREATE colData + DESeqDataSet####
# Build DESeq2 dataset and run the DESeq2 pipeline.

coldata <- samples %>%
  dplyr::select(sample, stage) %>%
  tibble::column_to_rownames("sample")  # DESeq2 requires rownames = sample IDs

dds <- DESeqDataSetFromTximport(
  txi,
  colData = coldata,
  design = ~ stage                      # model expression as a function of stage
)

# Filter out very low-count genes to reduce noise and improve estimation stability.
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2: size factors, dispersions, model fit, Wald tests.
dds <- DESeq(dds)

# Print coefficient names.
# With "early" as the reference, DESeq2 provides:
#   stage_thin_vs_early and stage_mature_vs_early
resultsNames(dds)


####7) OUTPUT FOLDERS: CREATE DIRECTORIES FOR RESULTS + FIGURES####
# Make sure output folders exist before writing files.

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)
dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)
dir.create("results/ora", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)


####8) FIGURE 1: PCA ON VST-TRANSFORMED COUNTS####
# VST improves interpretability by stabilizing variance across expression levels.
# PCA summarizes overall structure/outliers and whether samples cluster by stage.

vsd <- vst(dds, blind = FALSE)  # blind=FALSE respects the design during transformation

pca_df <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)  # get PCA coordinates
percentVar <- round(100 * attr(pca_df, "percentVar"))          # % variance explained

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = stage)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

ggsave("figures/fig1_pca.png", p_pca, width = 6, height = 5, dpi = 300)


####9) DIFFERENTIAL EXPRESSION: 3 CONTRASTS + VOLCANO PLOTS####
# Compute DE for 3 biologically meaningful comparisons.

contrasts <- list(
  thin_vs_early   = c("stage", "thin", "early"),
  mature_vs_early = c("stage", "mature", "early"),
  mature_vs_thin  = c("stage", "mature", "thin")
)

# res_list stores DE tables (one per contrast) so later chunks can reuse them.
res_list <- list()

# coef_map defines which contrasts correspond to existing model coefficients (shrinkable with apeglm).
# "mature_vs_thin" is NOT a model coefficient with early as reference, so it is computed directly via contrast.
coef_map <- c(
  thin_vs_early   = "stage_thin_vs_early",
  mature_vs_early = "stage_mature_vs_early"
)

for (nm in names(contrasts)) {
  
  # Pull the current contrast definition (factor name, numerator level, denominator level)
  con <- contrasts[[nm]]
  
  # If the comparison exists as a coefficient, use apeglm shrinkage for cleaner log2FC estimates.
  if (nm %in% names(coef_map)) {
    
    # lfcShrink reduces noisy log2FC estimates while generally preserving large effects.
    # coef=... must match one of resultsNames(dds).
    res_shrunk <- lfcShrink(dds, coef = coef_map[[nm]], type = "apeglm")
    
    # Convert DESeqResults object -> data.frame -> tibble; keep gene IDs; sort by padj.
    res_tbl <- as.data.frame(res_shrunk) %>%
      rownames_to_column("gene_id") %>%
      as_tibble() %>%
      arrange(padj)
    
  } else {
    
    # For contrasts not represented by a single coefficient (mature vs thin),
    # compute standard DESeq2 results using contrast=...
    res <- results(dds, contrast = con)
    
    res_tbl <- as.data.frame(res) %>%
      rownames_to_column("gene_id") %>%
      as_tibble() %>%
      arrange(padj)
  }
  
  # Save the DE table to disk for reporting/reproducibility.
  write_csv(res_tbl, file.path("results/deseq2", paste0("deseq2_", nm, ".csv")))
  
  # Save the tibble in memory for downstream steps (annotation, ORA, figure building).
  res_list[[nm]] <- res_tbl
  
  # Open PNG device for volcano plot output.
  png(file.path("figures", paste0("fig2_volcano_", nm, ".png")),
      width = 1800, height = 1400, res = 250)
  
  # Generate a volcano plot: x=log2FC, y=padj; thresholds p<0.05 and |log2FC|>=1.
  # Plot formatting is tuned to reduce wasted space from title/legend.
  p_volc <- EnhancedVolcano(
    res_tbl,
    lab = res_tbl$gene_id,        
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1.0,
    title = nm,
    labSize = 3,
    legendPosition = "right",
    legendLabSize = 10,
    legendIconSize = 3,
    titleLabSize = 14
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      plot.margin = margin(5, 5, 5, 5)   # tight outer margins so data area is larger
    )
  
  # Print the plot to the PNG device, then close device.
  print(p_volc)
  dev.off()
}


####10) FIGURE 3: HEATMAP OF TOP 50 MOST VARIABLE GENES (VST)####
# Select the 50 genes with the highest variance across samples (on VST scale)
# and plot a clustered heatmap to visualize global expression patterns.

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50) # indices of top variable genes
mat <- assay(vsd)[topVarGenes, , drop = FALSE]                        # VST expression matrix subset
mat <- mat - rowMeans(mat)                                            # mean-center each gene (row)

png("figures/fig3_heatmap_top50var.png", width = 1800, height = 1600, res = 250)
pheatmap(
  mat,
  annotation_col = as.data.frame(coldata["stage"]), # annotate columns by stage
  show_rownames = FALSE,                            # remove gene labels per instructor feedback
  show_colnames = TRUE
)
dev.off()


####11) FUNCTIONAL ANNOTATION (gene-level)####
# Functional annotation = mapping gene IDs to biological labels (names/descriptions/domains/GO/pathways).

keytype_use <- "ORF"  # DE gene IDs are yeast ORFs like YGR087C, so ORF is the correct key type.

annotate_genes <- function(gene_ids, keytype_use = "ORF") {
  # Purpose: given a vector of gene IDs, return a compact annotation table.
  # Retrieve:
  #   - basic fields (GENENAME/COMMON/DESCRIPTION/SGD/UNIPROT)
  #   - multi-valued fields (GO/ONTOLOGY/PATH/PFAM/INTERPRO) collapsed into semicolon lists
  
  gene_ids <- unique(gene_ids)  # avoid redundant lookups
  
  # ---- 1) Basic 1:1 annotation ----
  basic_cols <- c("GENENAME", "COMMON", "DESCRIPTION", "SGD", "UNIPROT")
  
  anno_basic <- AnnotationDbi::select(
    x       = org.Sc.sgd.db,   # annotation database
    keys    = gene_ids,        # keys to map
    keytype = keytype_use,     # what the keys represent (ORF IDs)
    columns = basic_cols       # which annotations to retrieve
  ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(gene_id = !!keytype_use) %>%   # unify key column name
    dplyr::group_by(gene_id) %>%                 # one output row per gene
    dplyr::summarise(
      # first() picks the first available non-missing value for each annotation column
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
  
  # ---- 4) Merge annotation tables ----
  # Left-joins keep all genes from the basic table even if GO/PFAM/etc. are missing.
  anno_basic %>%
    dplyr::left_join(anno_go, by = "gene_id") %>%
    dplyr::left_join(anno_other, by = "gene_id")
}

# Build a "master" annotation table for every gene that appears in any DE result table.
all_genes <- unique(unlist(lapply(res_list, function(x) x$gene_id)))
ann_all <- annotate_genes(all_genes, keytype_use = keytype_use)
readr::write_csv(ann_all, "results/annotation/functional_annotation_master.csv")

# For each contrast, annotate ONLY the significant DE genes and save a contrast-specific annotated table.
for (nm in names(res_list)) {
  
  res_tbl <- res_list[[nm]]  # DE results for this contrast
  
  # Define significant genes using the same thresholds used elsewhere (padj + effect size).
  sig_genes <- res_tbl %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  # If no significant genes, skip writing an annotated table.
  if (length(sig_genes) == 0) next
  
  # Retrieve functional annotations for significant genes.
  ann_sig <- annotate_genes(sig_genes, keytype_use = keytype_use)
  
  # Merge annotations onto the DE results table (adds name/description/GO/etc. columns).
  res_annotated <- res_tbl %>%
    dplyr::left_join(ann_sig, by = "gene_id") %>%
    dplyr::arrange(padj)
  
  # Save annotated DE table (useful for interpreting top DE genes in Results/Discussion).
  readr::write_csv(
    res_annotated,
    file.path("results/annotation", paste0("deseq2_", nm, "_annotated.csv"))
  )
}

####12) FIGURE 4: Functional annotation classification (Top GO BP terms by frequency)####
# This block was an earlier approach: count how many significant DE genes map to each GO BP term.
# Ultimately chose not to use this figure (saved as figX...).

if (!requireNamespace("GO.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("GO.db")
}
library(GO.db)

dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

keytype_use <- "ORF"

get_go_bp_counts <- function(sig_genes, contrast_name, keytype_use = "ORF") {
  # Purpose: for a set of significant genes, count how many unique genes map to each GO BP term.
  # Returns: a tibble with GO ID, n_genes, term name, and the contrast label for faceting.
  
  sig_genes <- unique(sig_genes)
  if (length(sig_genes) == 0) return(tibble())
  
  # Retrieve gene -> GO mappings (GO + ONTOLOGY together keeps category info aligned).
  go_map <- AnnotationDbi::select(
    x       = org.Sc.sgd.db,
    keys    = sig_genes,
    keytype = keytype_use,
    columns = c("GO", "ONTOLOGY")
  ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(gene_id = !!keytype_use) %>%
    dplyr::filter(!is.na(GO), !is.na(ONTOLOGY)) %>%   # drop missing mappings
    dplyr::filter(ONTOLOGY == "BP") %>%               # keep only Biological Process
    dplyr::distinct(gene_id, GO, .keep_all = TRUE)    # each gene counted once per GO term
  
  if (nrow(go_map) == 0) return(tibble())
  
  # Count number of genes per GO term (frequency among significant DE genes).
  counts <- go_map %>%
    dplyr::count(GO, name = "n_genes") %>%
    dplyr::arrange(dplyr::desc(n_genes))
  
  # Convert GO IDs to readable GO term names using GO.db.
  term_map <- AnnotationDbi::select(
    x       = GO.db,
    keys    = counts$GO,
    keytype = "GOID",
    columns = c("TERM")
  ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(GO = GOID)
  
  counts %>%
    dplyr::left_join(term_map, by = "GO") %>%
    dplyr::mutate(
      contrast = contrast_name,
      term_label = ifelse(is.na(TERM), GO, TERM)  # fallback to GO ID if term name missing
    )
}

# Build the GO BP frequency table for each contrast and combine into one data frame.
go_bp_all <- purrr::map_dfr(names(res_list), function(nm) {
  
  res_tbl <- res_list[[nm]]
  
  sig_genes <- res_tbl %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  get_go_bp_counts(sig_genes, contrast_name = nm, keytype_use = keytype_use)
})

readr::write_csv(go_bp_all, "results/annotation/go_bp_functional_classification_counts.csv")

# Select top N terms per contrast (highest n_genes) for plotting.
top_n <- 15
go_bp_top <- go_bp_all %>%
  dplyr::group_by(contrast) %>%
  dplyr::slice_max(order_by = n_genes, n = top_n, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(term_label = forcats::fct_reorder(term_label, n_genes)) %>% # order bars nicely
  dplyr::ungroup()

# Plot only if there is data available.
if (nrow(go_bp_top) == 0) {
  warning("No GO BP annotations found for significant genes under your thresholds; fig5 not generated.")
} else {
  
  # Horizontal bar chart of GO BP term frequencies by contrast.
  p_go_class <- ggplot2::ggplot(go_bp_top, ggplot2::aes(x = term_label, y = n_genes)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ contrast, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 11, margin = ggplot2::margin(t = 2, r = 8, b = 2, l = 8)),
      plot.margin = ggplot2::margin(t = 10, r = 25, b = 10, l = 10),
      plot.title = ggplot2::element_text(size = 16, hjust = 0.5)
    ) +
    ggplot2::labs(
      title = "Functional classification of significant DE genes (Top GO Biological Process terms)",
      x = "GO Biological Process term",
      y = "Number of significant DE genes annotated to term"
    )
  
  ggplot2::ggsave(
    filename = "figures/figX_functional_classification_topGO_BP.png",
    plot = p_go_class,
    width = 14,
    height = 7,
    dpi = 300,
    units = "in",
    limitsize = FALSE
  )
}

####13) ORA: GO BIOLOGICAL PROCESS ENRICHMENT (clusterProfiler)####
# ORA tests whether GO BP terms are over-represented among significant DE genes
# compared to the background set (all genes in the OrgDb / tested set).

for (nm in names(res_list)) {
  
  res_tbl <- res_list[[nm]]  # DE results for contrast nm
  
  # Define significant genes for ORA input (same thresholds used in Results).
  sig <- res_tbl %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  # If too few genes, enrichment will be unstable/uninformative, so skip.
  if (length(sig) < 10) next
  
  # Run enrichGO for Biological Process terms.
  ego <- enrichGO(
    gene = sig,
    OrgDb = org.Sc.sgd.db,
    keyType = keytype_use,     # ORF IDs must match the DE gene_id format
    ont = "BP",
    pAdjustMethod = "BH",      # Benjamini-Hochberg correction for enriched terms
    readable = FALSE           # avoids SYMBOL conversion issues in this OrgDb setup
  )
  
  # If no enriched terms returned, skip writing/plotting.
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next
  
  # Save ORA results table for this contrast.
  write_csv(as.data.frame(ego), file.path("results/ora", paste0("ora_go_bp_", nm, ".csv")))
  
  # Dotplot visualizes enriched terms with gene ratios and adjusted p-values.
  p <- dotplot(ego, showCategory = 15) + ggtitle(paste0("GO BP ORA: ", nm))
  
  # Save ORA dotplot figure.
  ggsave(
    file.path("figures", paste0("fig4_ora_go_bp_", nm, ".png")),
    p,
    width = 8, height = 6, dpi = 300
  )
}

####14) FIGURE 5: Functional classification of significant DE genes (PERCENTAGES)####
# This figure is a “functional classification” summary built from ORA output:
#   - Take the top ORA terms (lowest adjusted p-values) per contrast
#   - For each term, Count = number of significant genes in that GO term
#   - Convert Count to % of all significant genes in that contrast (denominator = sig_total)
#   - Remove the generic "biological_process" term if present

ora_files <- list.files("results/ora", pattern = "^ora_go_bp_.*\\.csv$", full.names = TRUE)
stopifnot(length(ora_files) > 0)  # ensures ORA results exist before proceeding

# Read ORA CSVs into a named list keyed by contrast name (e.g., "thin_vs_early")
ora_list <- ora_files %>%
  setNames(gsub("^ora_go_bp_|\\.csv$", "", basename(.))) %>%
  lapply(\(f) readr::read_csv(f, show_col_types = FALSE))

top_n_terms <- 10  # number of most significant ORA terms per contrast to display

# Build a table of the top ORA terms per contrast (used as the "term whitelist" for plotting).
top_terms_df <- purrr::imap_dfr(ora_list, function(df, nm) {
  df %>%
    dplyr::arrange(p.adjust) %>%            # rank by adjusted p-value
    dplyr::slice_head(n = top_n_terms) %>%  # take top N enriched terms
    dplyr::transmute(
      contrast = nm,
      go_id = ID,
      term = Description
    )
})

# Remove uninformative umbrella term (can appear depending on GO structure).
top_terms_df <- top_terms_df %>%
  dplyr::filter(!is.na(term), term != "biological_process")

# Build plotting data: one row per (contrast, term) with counts and percent of significant genes.
class_df <- purrr::imap_dfr(ora_list, function(df, nm) {
  
  # Denominator: total number of significant DE genes in this contrast.
  sig_total <- res_list[[nm]] %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    dplyr::pull(gene_id) %>%
    unique() %>%
    length()
  
  if (is.na(sig_total) || sig_total == 0) return(tibble::tibble())
  
  # Only keep the terms selected as “top terms” for this contrast.
  top_terms_nm <- top_terms_df %>% dplyr::filter(contrast == nm)
  
  # Join ORA results to the chosen terms and compute percent.
  df %>%
    dplyr::inner_join(
      top_terms_nm,
      by = c("ID" = "go_id", "Description" = "term")
    ) %>%
    dplyr::mutate(
      n_genes   = Count,                    # clusterProfiler Count column = DE genes in term
      percent   = 100 * (n_genes / sig_total),
      sig_total = sig_total,
      contrast  = nm,
      term      = Description
    ) %>%
    dplyr::filter(term != "biological_process") %>%
    dplyr::select(contrast, term, n_genes, sig_total, percent)
})

# Plot: horizontal bars of % significant genes per top ORA term, faceted by contrast.
p_fig5 <- ggplot(class_df, aes(x = percent, y = reorder(term, percent))) +
  geom_col() +
  facet_wrap(~ contrast, scales = "free_y") +
  labs(
    title = "Functional classification of significant DE genes (Top GO BP terms)",
    x = "Percent of significant DE genes annotated to term",
    y = "GO term"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 9)
  )

ggsave("figures/fig5_functional_classification_topGO_BP_percent.png",
       p_fig5, width = 12, height = 6, dpi = 300)

####15) FIGURE 6: 3-GENE EXPRESSION ACROSS EARLY, MEDIUM, AND LATE TIME POINTS####
# This figure visualizes normalized expression across stages for 3 genes of interest.
# Each stage shows 3 points because there are 3 biological replicates per stage.

genes_of_interest <- c("YIR019C", "YJL159W", "YJR151C")  # FLO11, HSP150, DAN4

# Sanity-check: pull log2FC and padj for these genes across all contrasts.
# This confirms that chosen genes are actually changing and supports interpretation in the Results.
check_goi <- purrr::imap_dfr(res_list, function(tbl, nm) {
  tbl %>%
    dplyr::filter(gene_id %in% genes_of_interest) %>%
    dplyr::mutate(contrast = nm) %>%
    dplyr::select(contrast, gene_id, log2FoldChange, padj) %>%
    dplyr::arrange(contrast, padj)
})
check_goi

# Build plotting data by extracting per-sample normalized counts for each gene.
# plotCounts(..., returnData=TRUE) returns a data.frame with:
#   - count: normalized count (DESeq2)
#   - stage: group label
#   - sample: sample ID (and other columns depending on DESeq2 version)
plot_df <- purrr::map_dfr(genes_of_interest, function(g) {
  d <- plotCounts(dds, gene = g, intgroup = "stage", returnData = TRUE)
  d$gene_id <- g  # tag rows so we can facet by gene later
  d
})

# Plot: points = each replicate, jittered to avoid overlap; line/points = mean per stage.
p <- ggplot(plot_df, aes(x = stage, y = count)) +
  geom_point(position = position_jitter(width = 0.12), size = 2) +  # replicate-level variation
  stat_summary(fun = mean, geom = "line", aes(group = 1)) +         # mean trend across stages
  stat_summary(fun = mean, geom = "point", size = 3) +             # mean marker per stage
  facet_wrap(~ gene_id, scales = "free_y") +                        # separate panel per gene
  theme_bw() +
  xlab("Stage") +
  ylab("Normalized expression (DESeq2)")

ggsave("figures/fig6_three_genes.png", p, width = 8, height = 4.5, dpi = 300)

####16) REPRODUCIBILITY: SAVE SESSION INFO####

writeLines(capture.output(sessionInfo()), "results/sessionInfo_R.txt")