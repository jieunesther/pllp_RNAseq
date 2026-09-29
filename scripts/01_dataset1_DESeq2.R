# ============================================================================
# 01_dataset1_DESeq2.R
# Dataset 1 (GEO GSE320085, 2021-03 run): EK / pd1116 / pd1245  x  HP / LP,
# 18 samples, HTSeq per-gene counts (Ensembl GRCz11, release 100).
#
# Produces everything downstream scripts need:
#   results/dataset1/raw_counts.csv               32,520 genes x 18 (Ensembl IDs)
#   results/dataset1/gene_annotation.csv          Ensembl ID -> symbol used as row name
#   results/dataset1/sample_metadata.csv
#   results/dataset1/normalized_counts.csv        DESeq2 median-of-ratios (22,863 genes)
#   results/dataset1/vst_transformed_counts.csv   VST (blind = FALSE, full model)
#   results/dataset1/PCA_coordinates.csv          plotPCA, top 500 variable genes (Fig 2B)
#   results/dataset1/DESeq2_results_<contrast>.csv       all genes, 9 contrasts
#   results/dataset1/DESeq2_significant_<contrast>.csv   padj < 0.05 & |LFC| > 1  (S1 Table)
#   results/dataset1/interaction_<allele>_genotype_x_diet.csv  (+ _significant)
#   results/dataset1/DEG_summary_table.csv        Fig 2C counts
#   results/dataset1/genes_unique_to_<allele>_<diet>_table.csv  (S2 Table)
#   results/dataset1/DEG_overlaps.csv             Venn counts (Fig S2A-C)
#   results/dataset1/QC_*.pdf, PCA_all_samples.pdf, DEG_summary_plot.pdf, venn_diagrams.pdf
#
# Design (Methods): two models on the same filtered counts
#   full  : ~ genotype + diet + genotype:diet   -> VST (PCA, heatmaps), interaction terms
#   group : ~ condition (genotype_diet)          -> all pairwise Wald contrasts
# Filter : genes with >= 10 counts in >= 3 samples (32,520 -> 22,863)
# DEG    : BH-adjusted p < 0.05 and |log2FC| > 1
# ============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), "..", "R", "common.R"))
need_pkgs("DESeq2", "org.Dr.eg.db", "AnnotationDbi", "ggplot2", "VennDiagram", "grid")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
OUT <- RES_D1

# ============================================================================
# 1. SAMPLE METADATA
# ============================================================================
banner("1. Sample metadata")
sample_info <- read.csv(file.path(META, "dataset1_samples.csv"), stringsAsFactors = FALSE)
sample_info$genotype  <- factor(sample_info$genotype, levels = c("EK", "pd1116", "pd1245"))
sample_info$diet      <- factor(sample_info$diet, levels = c("HP", "LP"))
sample_info$condition <- factor(paste(sample_info$genotype, sample_info$diet, sep = "_"), levels = COND_LEVELS)
print(table(sample_info$genotype, sample_info$diet))
stopifnot(nrow(sample_info) == 18, all(table(sample_info$condition) == 3))
write.csv(sample_info, file.path(OUT, "sample_metadata.csv"), row.names = FALSE)

# ============================================================================
# 2. IMPORT HTSeq COUNTS
# ============================================================================
# GEO supplies one HTSeq count file per sample. Files are matched to samples by
# sample_id appearing in the file name (GEO prefixes files with GSM IDs), and
# may be gzipped. Alternatively a single matrix data/dataset1_htseq_counts/
# raw_counts.csv (gene_id column + one column per sample_id) is accepted.
banner("2. Import HTSeq counts")
count_dir <- file.path(DATA, "dataset1_htseq_counts")
files <- list.files(count_dir, full.names = TRUE)
matrix_file <- files[basename(files) %in% c("raw_counts.csv", "raw_counts.csv.gz")]

read_count_file <- function(f, id) {
  con <- if (grepl("\\.gz$", f)) gzfile(f) else f
  d <- read.table(con, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("gene_id", id), quote = "", comment.char = "")
  d
}

if (length(matrix_file)) {
  kv("input", basename(matrix_file[1]))
  count_data <- read.csv(matrix_file[1], row.names = 1, check.names = FALSE)
} else {
  count_data <- NULL
  for (i in seq_len(nrow(sample_info))) {
    id <- sample_info$sample_id[i]; gsm <- sample_info$gsm[i]
    # match by GEO accession (GSM9533663_...) or by sample title (EK_HP_1_S1)
    f <- files[!grepl("README", basename(files)) &
               (grepl(gsm, basename(files), fixed = TRUE) | grepl(id, basename(files), fixed = TRUE))]
    if (length(f) != 1) stop("Expected exactly one count file for ", id, " (", gsm, ") in ", count_dir,
                             " (found ", length(f), "). See data/dataset1_htseq_counts/README.md")
    kv(id, basename(f))
    d <- read_count_file(f, id)
    count_data <- if (is.null(count_data)) d else merge(count_data, d, by = "gene_id")
  }
  count_data <- count_data[!grepl("^__", count_data$gene_id), ]      # drop HTSeq summary rows
  rownames(count_data) <- count_data$gene_id
  count_data <- count_data[, -1]
}
count_data <- as.matrix(count_data[, sample_info$sample_id])
storage.mode(count_data) <- "integer"
kv("raw count matrix", paste(dim(count_data), collapse = " x "))
kv("library sizes (M reads)", paste(round(colSums(count_data) / 1e6, 1), collapse = " "))
check_expected("d1_raw_genes", nrow(count_data))
check_expected("d1_samples", ncol(count_data))
write.csv(count_data, file.path(OUT, "raw_counts.csv"), quote = FALSE)

pdf(file.path(OUT, "QC_library_sizes.pdf"), width = 10, height = 6)
print(ggplot(data.frame(sample_info, total = colSums(count_data)),
             aes(x = sample_id, y = total / 1e6, fill = condition)) +
        geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Sample", y = "Assigned reads (millions)", title = "Library sizes"))
invisible(dev.off())

# ============================================================================
# 3. LOW-COUNT FILTER AND GENE SYMBOLS
# ============================================================================
banner("3. Filter (>=10 counts in >=3 samples) and symbol annotation")
keep <- rowSums(count_data >= 10) >= 3
kv("genes before / after filter", c(nrow(count_data), sum(keep)))
check_expected("d1_filtered_genes", sum(keep))
counts_f <- count_data[keep, ]

# Ensembl -> symbol via org.Dr.eg.db (first hit). Unmapped genes keep their
# Ensembl ID; duplicated symbols are made unique (.1, .2 ...). The mapping used
# is written out so every downstream table can be traced back to an Ensembl ID
# regardless of the annotation-package version.
sym <- AnnotationDbi::mapIds(org.Dr.eg.db, keys = rownames(counts_f), column = "SYMBOL",
                             keytype = "ENSEMBL", multiVals = "first")
sym <- ifelse(is.na(sym), rownames(counts_f), sym)
sym_used <- make.unique(as.character(sym))
annot <- data.frame(ensembl_id = rownames(counts_f), symbol = sym, symbol_used = sym_used,
                    stringsAsFactors = FALSE)
kv("genes with a symbol", sum(annot$symbol != annot$ensembl_id))
kv("genes kept as Ensembl ID", sum(annot$symbol == annot$ensembl_id))
kv("symbols de-duplicated", sum(annot$symbol != annot$symbol_used))
kv("org.Dr.eg.db version", as.character(packageVersion("org.Dr.eg.db")))
write.csv(annot, file.path(OUT, "gene_annotation.csv"), row.names = FALSE)
rownames(counts_f) <- sym_used

# ============================================================================
# 4. FULL MODEL: normalisation, VST, interaction terms
# ============================================================================
banner("4. DESeq2 full model  ~ genotype + diet + genotype:diet")
dds <- DESeqDataSetFromMatrix(counts_f, colData = sample_info, design = ~ genotype + diet + genotype:diet)
dds <- DESeq(dds, quiet = TRUE)
kv("resultsNames", resultsNames(dds))

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(OUT, "normalized_counts.csv"), quote = FALSE)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)
write.csv(vsd_mat, file.path(OUT, "vst_transformed_counts.csv"), quote = FALSE)
kv("VST matrix", paste(dim(vsd_mat), collapse = " x "))
kv("VST spot-check [1,1]", paste(rownames(vsd_mat)[1], round(vsd_mat[1, 1], 4)))

# Genotype x diet interaction: is the LP-vs-HP effect different in the mutant than in EK?
for (allele in c("pd1116", "pd1245")) {
  term <- paste0("genotype", allele, ".dietLP")
  res  <- results(dds, name = term, alpha = ALPHA)
  df   <- as.data.frame(res[order(res$padj), ]); df$gene_id <- rownames(df)
  df$ensembl_id <- annot$ensembl_id[match(df$gene_id, annot$symbol_used)]
  write.csv(df, file.path(OUT, paste0("interaction_", sub("pd", "", allele), "_genotype_x_diet.csv")), row.names = FALSE)
  sig <- df[!is.na(df$padj) & df$padj < ALPHA & abs(df$log2FoldChange) > LFC, ]
  write.csv(sig, file.path(OUT, paste0("interaction_", sub("pd", "", allele), "_significant.csv")), row.names = FALSE)
  kv(paste("interaction", term, "DEGs"), nrow(sig))
}

# ============================================================================
# 5. PCA (Fig 2B): plotPCA on the 500 most variable VST genes
# ============================================================================
banner("5. PCA")
pcaData <- plotPCA(vsd, intgroup = c("genotype", "diet"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
kv("PC1 / PC2 % variance", percentVar[1:2])
check_expected("PCA_PC1_pct", percentVar[1]); check_expected("PCA_PC2_pct", percentVar[2])
write.csv(pcaData, file.path(OUT, "PCA_coordinates.csv"), row.names = FALSE)
pdf(file.path(OUT, "PCA_all_samples.pdf"), width = 7, height = 6)
print(ggplot(pcaData, aes(PC1, PC2, color = genotype, shape = diet)) + geom_point(size = 4) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        coord_fixed() + theme_bw() + scale_color_manual(values = GENO_COLORS))
invisible(dev.off())

# ============================================================================
# 6. GROUP MODEL: nine pairwise contrasts (Fig 2C, S1 Table)
# ============================================================================
banner("6. DESeq2 group model  ~ condition, 9 Wald contrasts")
dds_cond <- DESeqDataSetFromMatrix(counts_f, colData = sample_info, design = ~ condition)
dds_cond <- DESeq(dds_cond, quiet = TRUE)

comparisons <- list(
  HP_1116_vs_EK   = c("pd1116_HP", "EK_HP"),
  HP_1245_vs_EK   = c("pd1245_HP", "EK_HP"),
  HP_1245_vs_1116 = c("pd1245_HP", "pd1116_HP"),
  LP_1116_vs_EK   = c("pd1116_LP", "EK_LP"),
  LP_1245_vs_EK   = c("pd1245_LP", "EK_LP"),
  LP_1245_vs_1116 = c("pd1245_LP", "pd1116_LP"),
  EK_HP_vs_LP     = c("EK_HP", "EK_LP"),
  `1116_HP_vs_LP` = c("pd1116_HP", "pd1116_LP"),
  `1245_HP_vs_LP` = c("pd1245_HP", "pd1245_LP"))
labels <- c("HP: pd1116 vs EK", "HP: pd1245 vs EK", "HP: pd1245 vs pd1116",
            "LP: pd1116 vs EK", "LP: pd1245 vs EK", "LP: pd1245 vs pd1116",
            "EK: HP vs LP", "pd1116: HP vs LP", "pd1245: HP vs LP")

all_comparisons <- list(); summary_rows <- list()
for (nm in names(comparisons)) {
  cc  <- comparisons[[nm]]
  res <- results(dds_cond, contrast = c("condition", cc[1], cc[2]), alpha = ALPHA)
  df  <- as.data.frame(res[order(res$padj), ]); df$gene_id <- rownames(df)
  df$ensembl_id <- annot$ensembl_id[match(df$gene_id, annot$symbol_used)]
  write.csv(df, file.path(OUT, paste0("DESeq2_results_", nm, ".csv")), row.names = FALSE)
  sig <- df[!is.na(df$padj) & df$padj < ALPHA & abs(df$log2FoldChange) > LFC, ]
  write.csv(sig, file.path(OUT, paste0("DESeq2_significant_", nm, ".csv")), row.names = FALSE)
  up <- sum(sig$log2FoldChange > 0); dn <- sum(sig$log2FoldChange < 0)
  cat(sprintf("  %-16s %s vs %s : %4d DEGs (%3d up / %3d down)\n", nm, cc[1], cc[2], nrow(sig), up, dn))
  check_expected(paste0("DEG_", nm), nrow(sig))
  check_expected(paste0("DEG_", nm, "_up"), up)
  check_expected(paste0("DEG_", nm, "_down"), dn)
  all_comparisons[[nm]] <- list(all = df, sig = sig)
  summary_rows[[nm]] <- data.frame(Comparison = labels[match(nm, names(comparisons))],
                                   Total_DEGs = nrow(sig), Upregulated = up, Downregulated = dn)
}
summary_table <- do.call(rbind, summary_rows); rownames(summary_table) <- NULL
write.csv(summary_table, file.path(OUT, "DEG_summary_table.csv"), row.names = FALSE)
print(summary_table)

deg_union <- unique(unlist(lapply(all_comparisons, function(x) x$sig$gene_id)))
kv("union of DEGs across 9 contrasts (Fig 2A)", length(deg_union))
check_expected("DEG_union_9_comparisons", length(deg_union))
writeLines(deg_union, file.path(OUT, "DEG_union_9_comparisons.txt"))

# ============================================================================
# 7. OVERLAPS (Fig S2A-C, S2 Table)
# ============================================================================
banner("7. DEG overlaps")
S <- lapply(all_comparisons, function(x) x$sig$gene_id)
ov <- c(
  venn_diet_all3         = length(Reduce(intersect, S[c("EK_HP_vs_LP", "1116_HP_vs_LP", "1245_HP_vs_LP")])),
  venn_diet_EK_and_1116  = length(intersect(S$EK_HP_vs_LP, S$`1116_HP_vs_LP`)),
  venn_diet_EK_and_1245  = length(intersect(S$EK_HP_vs_LP, S$`1245_HP_vs_LP`)),
  venn_diet_unique_EK    = length(setdiff(S$EK_HP_vs_LP, union(S$`1116_HP_vs_LP`, S$`1245_HP_vs_LP`))),
  venn_diet_unique_1116  = length(setdiff(S$`1116_HP_vs_LP`, union(S$EK_HP_vs_LP, S$`1245_HP_vs_LP`))),
  venn_diet_unique_1245  = length(setdiff(S$`1245_HP_vs_LP`, union(S$EK_HP_vs_LP, S$`1116_HP_vs_LP`))),
  venn_HP_shared         = length(intersect(S$HP_1116_vs_EK, S$HP_1245_vs_EK)),
  venn_HP_unique_1116    = length(setdiff(S$HP_1116_vs_EK, S$HP_1245_vs_EK)),
  venn_HP_unique_1245    = length(setdiff(S$HP_1245_vs_EK, S$HP_1116_vs_EK)),
  venn_LP_shared         = length(intersect(S$LP_1116_vs_EK, S$LP_1245_vs_EK)),
  venn_LP_unique_1116    = length(setdiff(S$LP_1116_vs_EK, S$LP_1245_vs_EK)),
  venn_LP_unique_1245    = length(setdiff(S$LP_1245_vs_EK, S$LP_1116_vs_EK)))
for (k in names(ov)) check_expected(k, ov[[k]])
write.csv(data.frame(overlap = names(ov), n = unname(ov)), file.path(OUT, "DEG_overlaps.csv"), row.names = FALSE)

# Genotype-specific DEG tables (S2 Table): DEGs vs EK in one allele but not the other
uniq <- list(
  genes_unique_to_1116_HP_table = all_comparisons$HP_1116_vs_EK$sig[!all_comparisons$HP_1116_vs_EK$sig$gene_id %in% S$HP_1245_vs_EK, ],
  genes_unique_to_1116_LP_table = all_comparisons$LP_1116_vs_EK$sig[!all_comparisons$LP_1116_vs_EK$sig$gene_id %in% S$LP_1245_vs_EK, ],
  genes_unique_to_1245_HP_table = all_comparisons$HP_1245_vs_EK$sig[!all_comparisons$HP_1245_vs_EK$sig$gene_id %in% S$HP_1116_vs_EK, ],
  genes_unique_to_1245_LP_table = all_comparisons$LP_1245_vs_EK$sig[!all_comparisons$LP_1245_vs_EK$sig$gene_id %in% S$LP_1116_vs_EK, ])
for (nm in names(uniq)) { write.csv(uniq[[nm]], file.path(OUT, paste0(nm, ".csv")), row.names = FALSE); kv(nm, nrow(uniq[[nm]])) }

# Quick-look Venn diagrams (publication versions are drawn in 08_figureS2.R)
pdf(file.path(OUT, "venn_diagrams.pdf"), width = 8, height = 8)
vd <- function(x, fills, main) { grid.newpage(); grid.draw(venn.diagram(x, filename = NULL, fill = fills, alpha = .5,
                                   scaled = FALSE, main = main, print.mode = c("raw", "percent"), margin = .1)) }
vd(list(pd1116 = S$HP_1116_vs_EK, pd1245 = S$HP_1245_vs_EK), GENO_COLORS[2:3], "HP: DEGs vs EK")
vd(list(pd1116 = S$LP_1116_vs_EK, pd1245 = S$LP_1245_vs_EK), GENO_COLORS[2:3], "LP: DEGs vs EK")
vd(list(EK = S$EK_HP_vs_LP, pd1116 = S$`1116_HP_vs_LP`, pd1245 = S$`1245_HP_vs_LP`), GENO_COLORS, "HP vs LP within genotype")
invisible(dev.off())

# Quick-look DEG summary bar chart (publication version: 05_figure2.R)
pd <- summary_table
pd$Category <- ifelse(grepl("HP vs LP", pd$Comparison), "Diet effects", "Genotype effects")
pd$Comparison <- factor(pd$Comparison, levels = rev(pd$Comparison))
pdf(file.path(OUT, "DEG_summary_plot.pdf"), width = 8, height = 6)
print(ggplot(pd) + geom_col(aes(Comparison, -Downregulated, fill = "Down"), width = .7) +
        geom_col(aes(Comparison, Upregulated, fill = "Up"), width = .7) +
        facet_grid(Category ~ ., scales = "free_y", space = "free_y") + coord_flip() +
        scale_y_continuous(labels = abs) + scale_fill_manual(values = c(Down = "cornflowerblue", Up = "firebrick3"), name = NULL) +
        theme_minimal() + labs(x = NULL, y = "Number of DEGs", title = "DEGs per comparison (FDR < 0.05, |LFC| > 1)"))
invisible(dev.off())

# ============================================================================
# 8. OUTPUT VERIFICATION
# ============================================================================
banner("8. Output verification")
verify_files(file.path(OUT, c("raw_counts.csv", "gene_annotation.csv", "sample_metadata.csv", "normalized_counts.csv",
                              "vst_transformed_counts.csv", "PCA_coordinates.csv", "DEG_summary_table.csv",
                              "DEG_overlaps.csv", "DEG_union_9_comparisons.txt",
                              paste0("DESeq2_significant_", names(comparisons), ".csv"),
                              paste0("DESeq2_results_", names(comparisons), ".csv"),
                              paste0(names(uniq), ".csv"),
                              "interaction_1116_genotype_x_diet.csv", "interaction_1245_genotype_x_diet.csv",
                              "QC_library_sizes.pdf", "PCA_all_samples.pdf", "DEG_summary_plot.pdf", "venn_diagrams.pdf")))
write_session()
cat("\nDone: 01_dataset1_DESeq2.R\n")
