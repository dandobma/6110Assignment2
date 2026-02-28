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
  
  library(tidyverse)
  library(tximport)
  library(DESeq2)
  library(matrixStats)
  library(pheatmap)
  library(EnhancedVolcano)
  
  library(clusterProfiler)
  library(org.Sc.sgd.db)
  library(AnnotationDbi)
  library(enrichplot)
})


####2) INPUTS: READ SAMPLE METADATA + DEFINE GROUPS####

samples <- readr::read_csv("data/samples.csv", show_col_types = FALSE) %>%
  mutate(stage = factor(stage, levels = c("early", "thin", "mature")))


####3) INPUTS: LOCATE SALMON QUANT FILES (quant.sf)#####

files <- file.path("results/salmon", samples$sample, "quant.sf")
names(files) <- samples$sample
stopifnot(all(file.exists(files)))


####4) BUILD tx2gene MAPPING FROM ENSEMBL GTF####

gtf <- "reference/Saccharomyces_cerevisiae.R64-1-1.111.gtf.gz"

gtf_tbl <- readr::read_tsv(
  gtf,
  comment = "#",
  col_names = FALSE,
  show_col_types = FALSE
)

colnames(gtf_tbl) <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

tx2gene <- gtf_tbl %>%
  dplyr::filter(feature == "transcript") %>%
  dplyr::transmute(
    transcript_id = stringr::str_match(attribute, 'transcript_id "([^"]+)"')[, 2],
    gene_id       = stringr::str_match(attribute, 'gene_id "([^"]+)"')[, 2]
  ) %>%
  dplyr::filter(!is.na(transcript_id), !is.na(gene_id)) %>%
  dplyr::distinct()


####5) tximport: IMPORT SALMON + SUMMARIZE TRANSCRIPT -> GENE####

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)


####6) DESeq2 SETUP: CREATE colData + DESeqDataSet####

coldata <- samples %>%
  dplyr::select(sample, stage) %>%
  tibble::column_to_rownames("sample")

dds <- DESeqDataSetFromTximport(
  txi,
  colData = coldata,
  design = ~ stage
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# Print coefficient names
resultsNames(dds)


####7) OUTPUT FOLDERS: CREATE DIRECTORIES FOR RESULTS + FIGURES####

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)
dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)
dir.create("results/ora", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)


####8) FIGURE 1: PCA ON VST-TRANSFORMED COUNTS####

vsd <- vst(dds, blind = FALSE)

pca_df <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = stage)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

ggsave("figures/fig1_pca.png", p_pca, width = 6, height = 5, dpi = 300)


####9) DIFFERENTIAL EXPRESSION: 3 CONTRASTS + VOLCANO PLOTS####

contrasts <- list(
  thin_vs_early   = c("stage", "thin", "early"),
  mature_vs_early = c("stage", "mature", "early"),
  mature_vs_thin  = c("stage", "mature", "thin")
)

# Store DE tables for later steps (annotation + ORA)
res_list <- list()

# Only coefficients available with early as reference:
coef_map <- c(
  thin_vs_early   = "stage_thin_vs_early",
  mature_vs_early = "stage_mature_vs_early"
)

for (nm in names(contrasts)) {
  
  con <- contrasts[[nm]]
  
  if (nm %in% names(coef_map)) {
    
    res_shrunk <- lfcShrink(dds, coef = coef_map[[nm]], type = "apeglm")
    
    res_tbl <- as.data.frame(res_shrunk) %>%
      rownames_to_column("gene_id") %>%
      as_tibble() %>%
      arrange(padj)
    
  } else {
    
    res <- results(dds, contrast = con)
    
    res_tbl <- as.data.frame(res) %>%
      rownames_to_column("gene_id") %>%
      as_tibble() %>%
      arrange(padj)
  }
  
  write_csv(res_tbl, file.path("results/deseq2", paste0("deseq2_", nm, ".csv")))
  res_list[[nm]] <- res_tbl
  
  png(file.path("figures", paste0("fig2_volcano_", nm, ".png")),
      width = 1800, height = 1400, res = 250)
  
#  print(EnhancedVolcano(
#    res_tbl,
#    lab = res_tbl$gene_id,
#    x = "log2FoldChange",
#    y = "padj",
#    pCutoff = 0.05,
#    FCcutoff = 1.0,
#    title = nm
#  ))
  
  p_volc <- EnhancedVolcano(
    res_tbl,
    lab = res_tbl$gene_id,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1.0,
    title = nm,
    labSize = 3,
    legendPosition = "right",     # move legend to side
    legendLabSize = 10,
    legendIconSize = 3,
    titleLabSize = 14              # reduce title size
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      plot.margin = margin(5, 5, 5, 5)   # tighten margins
    )
  
  print(p_volc)
  dev.off()
}


####10) FIGURE 3: HEATMAP OF TOP 50 MOST VARIABLE GENES (VST)####

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[topVarGenes, , drop = FALSE]
mat <- mat - rowMeans(mat)

png("figures/fig3_heatmap_top50var.png", width = 1800, height = 1600, res = 250)
pheatmap(mat, annotation_col = as.data.frame(coldata["stage"]), 
         show_rownames = F, show_colnames = T)
dev.off()


####11) FUNCTIONAL ANNOTATION (gene-level)####
#       Map ORF IDs -> gene names/descriptions + GO/PFAM/PATH annotations.

# DE gene IDs are yeast ORFs (e.g., YGR087C).
keytype_use <- "ORF"

annotate_genes <- function(gene_ids, keytype_use = "ORF") {
  
  gene_ids <- unique(gene_ids)
  
  # ---- 1) Basic 1:1 annotation ----
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
  
  anno_basic %>%
    dplyr::left_join(anno_go, by = "gene_id") %>%
    dplyr::left_join(anno_other, by = "gene_id")
}

# Save a master annotation table for all genes seen in DE results
all_genes <- unique(unlist(lapply(res_list, function(x) x$gene_id)))
ann_all <- annotate_genes(all_genes, keytype_use = keytype_use)
readr::write_csv(ann_all, "results/annotation/functional_annotation_master.csv")

# Annotate + save an annotated DE table for each contrast
for (nm in names(res_list)) {
  
  res_tbl <- res_list[[nm]]
  
  sig_genes <- res_tbl %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
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

####12) FIGURE 4: Functional annotation classification (Top GO BP terms by frequency)####
#        For each contrast:
#          - define significant DE genes (padj < 0.05, |log2FC| >= 1)
#          - pull GO annotations for those genes from org.Sc.sgd.db (GO + ONTOLOGY)
#          - keep BP only
#          - count unique genes per GO term
#          - map GO IDs -> human-readable term names using GO.db
#
#      Output:
#        figures/figX_functional_classification_topGO_BP.png (chose not to use)
#        results/annotation/go_bp_functional_classification_counts.csv

# GO term-name database (Bioconductor). Install once if needed.
if (!requireNamespace("GO.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("GO.db")
}
library(GO.db)

# Output folders
dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# Keytype must match gene IDs (YAL001C-like ORFs)
keytype_use <- "ORF"

# Helper: build GO BP term frequency table for one contrast
get_go_bp_counts <- function(sig_genes, contrast_name, keytype_use = "ORF") {
  sig_genes <- unique(sig_genes)
  if (length(sig_genes) == 0) return(tibble())
  
  # Pull gene -> GO mappings (keeps GO + ONTOLOGY paired correctly)
  go_map <- AnnotationDbi::select(
    x       = org.Sc.sgd.db,
    keys    = sig_genes,
    keytype = keytype_use,
    columns = c("GO", "ONTOLOGY")
  ) %>%
    tibble::as_tibble() %>%
    dplyr::rename(gene_id = !!keytype_use) %>%
    dplyr::filter(!is.na(GO), !is.na(ONTOLOGY)) %>%
    dplyr::filter(ONTOLOGY == "BP") %>%
    dplyr::distinct(gene_id, GO, .keep_all = TRUE)
  
  if (nrow(go_map) == 0) return(tibble())
  
  # Count how many unique genes map to each GO BP term
  counts <- go_map %>%
    dplyr::count(GO, name = "n_genes") %>%
    dplyr::arrange(dplyr::desc(n_genes))
  
  # Map GO IDs -> term names (TERM)
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
      term_label = ifelse(is.na(TERM), GO, TERM)
    )
}

# Build classification counts for all contrasts from res_list
# (res_list must exist from DE loop)
go_bp_all <- purrr::map_dfr(names(res_list), function(nm) {
  
  res_tbl <- res_list[[nm]]
  
  sig_genes <- res_tbl %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  get_go_bp_counts(sig_genes, contrast_name = nm, keytype_use = keytype_use)
})

# Save the full table
readr::write_csv(go_bp_all, "results/annotation/go_bp_functional_classification_counts.csv")

# Choose top N terms per contrast for plotting
top_n <- 15
go_bp_top <- go_bp_all %>%
  dplyr::group_by(contrast) %>%
  dplyr::slice_max(order_by = n_genes, n = top_n, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  # reorder terms within each facet for clean plotting
  dplyr::group_by(contrast) %>%
  dplyr::mutate(term_label = forcats::fct_reorder(term_label, n_genes)) %>%
  dplyr::ungroup()

# If nothing to plot, stop with a helpful message
if (nrow(go_bp_top) == 0) {
  warning("No GO BP annotations found for significant genes under your thresholds; fig5 not generated.")
} else {
  
  p_go_class <- ggplot2::ggplot(go_bp_top, ggplot2::aes(x = term_label, y = n_genes)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ contrast, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      # make facet strip text smaller
      strip.text = ggplot2::element_text(size = 11, margin = ggplot2::margin(t = 2, r = 8, b = 2, l = 8)),
      # a bit more spacing around the whole plot so nothing is clipped
      plot.margin = ggplot2::margin(t = 10, r = 25, b = 10, l = 10),
      # slightly reduce title size to avoid crowding
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
    width = 14,      # wider so facet strip labels fit
    height = 7,
    dpi = 300,
    units = "in",
    limitsize = FALSE
  )
}

####13) ORA: GO BIOLOGICAL PROCESS ENRICHMENT (clusterProfiler)####

for (nm in names(res_list)) {
  
  res_tbl <- res_list[[nm]]
  
  sig <- res_tbl %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  if (length(sig) < 10) next
  
  ego <- enrichGO(
    gene = sig,
    OrgDb = org.Sc.sgd.db,
    keyType = keytype_use,   # <-- FIXED: use ORF to match your gene_id column
    ont = "BP",
    pAdjustMethod = "BH",
    readable = FALSE         # <-- avoids SYMBOL conversion error for this OrgDb
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next
  
  write_csv(as.data.frame(ego), file.path("results/ora", paste0("ora_go_bp_", nm, ".csv")))
  
  p <- dotplot(ego, showCategory = 15) + ggtitle(paste0("GO BP ORA: ", nm))
  
  ggsave(
    file.path("figures", paste0("fig4_ora_go_bp_", nm, ".png")),
    p,
    width = 8, height = 6, dpi = 300
  )
}

####14) FIGURE 5: Functional classification of significant DE genes (PERCENTAGES)####
#       - For each contrast, classify significant DE genes into top GO BP terms
#       - Plot % of significant DE genes annotated to each term (within contrast)
#       - Remove the useless "biological_process" term

# Safely read ORA files if they exist
ora_files <- list.files("results/ora", pattern = "^ora_go_bp_.*\\.csv$", full.names = TRUE)
stopifnot(length(ora_files) > 0)

# Read ORA results into a named list by contrast (must match res_list names)
ora_list <- ora_files %>%
  setNames(gsub("^ora_go_bp_|\\.csv$", "", basename(.))) %>%
  lapply(\(f) readr::read_csv(f, show_col_types = FALSE))

# Choose how many GO BP terms to display per contrast
top_n_terms <- 10

# Build a data frame giving, per contrast, the top ORA terms to use for classification
top_terms_df <- purrr::imap_dfr(ora_list, function(df, nm) {
  df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n_terms) %>%
    dplyr::transmute(
      contrast = nm,
      go_id = ID,
      term = Description
    )
})

# Remove the useless generic category if it appears
top_terms_df <- top_terms_df %>%
  dplyr::filter(!is.na(term), term != "biological_process")

# ---- Figure 5 data: Functional classification as % of significant DE genes ----
class_df <- purrr::imap_dfr(ora_list, function(df, nm) {
  
  # denominator = number of significant DE genes in this contrast
  sig_total <- res_list[[nm]] %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    dplyr::pull(gene_id) %>%
    unique() %>%
    length()
  
  if (is.na(sig_total) || sig_total == 0) return(tibble::tibble())
  
  # term list for this contrast
  top_terms_nm <- top_terms_df %>% dplyr::filter(contrast == nm)
  
  # Join ORA results to the selected term list
  df %>%
    dplyr::inner_join(
      top_terms_nm,
      by = c("ID" = "go_id", "Description" = "term")
    ) %>%
    dplyr::mutate(
      n_genes  = Count,
      percent  = 100 * (n_genes / sig_total),
      sig_total = sig_total,
      contrast = nm,
      term = Description
    ) %>%
    dplyr::filter(term != "biological_process") %>%
    dplyr::select(contrast, term, n_genes, sig_total, percent)
})

# ---- Figure 5 plot ----
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


# Choose genes of interest
genes_of_interest <- c("YIR019C", "YJL159W", "YJR151C")  # FLO11, HSP150, DAN4

# check significance / effect size in each contrast
check_goi <- purrr::imap_dfr(res_list, function(tbl, nm) {
  tbl %>%
    dplyr::filter(gene_id %in% genes_of_interest) %>%
    dplyr::mutate(contrast = nm) %>%
    dplyr::select(contrast, gene_id, log2FoldChange, padj) %>%
    dplyr::arrange(contrast, padj)
})

check_goi

# pick “backup” candidates: top by padj among strong effects
#backup <- res_list[["mature_vs_early"]] %>%
#  dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
#  dplyr::arrange(padj) %>%
#  dplyr::slice_head(n = 30) %>%
#  dplyr::pull(gene_id)

#backup

plot_df <- purrr::map_dfr(genes_of_interest, function(g) {
  d <- plotCounts(dds, gene = g, intgroup = "stage", returnData = TRUE)
  d$gene_id <- g
  d
})

p <- ggplot(plot_df, aes(x = stage, y = count)) +
  geom_point(position = position_jitter(width = 0.12), size = 2) +
  stat_summary(fun = mean, geom = "line", aes(group = 1)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  facet_wrap(~ gene_id, scales = "free_y") +
  theme_bw() +
  xlab("Stage") +
  ylab("Normalized expression (DESeq2)")

ggsave("figures/fig6_three_genes.png", p, width = 8, height = 4.5, dpi = 300)

####16) REPRODUCIBILITY: SAVE SESSION INFO####


writeLines(capture.output(sessionInfo()), "results/sessionInfo_R.txt")
