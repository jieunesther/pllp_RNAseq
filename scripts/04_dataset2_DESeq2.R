# ============================================================================
# 04_dataset2_DESeq2.R
# Dataset 2 (GEO GSE320085, "TLHYNJ" run, Fig S7): EK vs adapted F8 pd1245,
# HP / LP, 14 samples. Input is the Plasmidsaurus featureCounts matrix
# (TLHYNJ-expression-matrix.tsv); only the *_count columns are used and
# fractional multi-mapping counts are rounded to integers for DESeq2.
#
# Same pipeline as dataset 1 (Methods): filter >=10 counts in >=3 samples;
#   full model  ~ genotype + diet + genotype:diet  -> VST (PCA Fig S7D, heatmap S7E)
#   group model ~ condition                         -> pd1245 vs EK within each diet
# DEG: BH-adjusted p < 0.05 and |log2FC| > 1.
#
# Outputs (results/dataset2/):
#   raw_counts.csv, gene_annotation.csv, sample_metadata.csv
#   vst_transformed_counts.csv, PCA_coordinates.csv, PCA_S7D.pdf
#   DESeq2_results_<diet>_pd1245_vs_EK.csv, DESeq2_significant_<diet>_pd1245_vs_EK.csv
# ============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), "..", "R", "common.R"))
need_pkgs("DESeq2", "ggplot2")
OUT <- RES_D2

banner("1. Input matrix")
mat_file <- list.files(file.path(DATA, "dataset2_TLHYNJ"), pattern = "expression-matrix\\.tsv(\\.gz)?$", full.names = TRUE)
if (length(mat_file) != 1) stop("Expected one *expression-matrix.tsv in data/dataset2_TLHYNJ/ (see its README.md)")
m <- read.delim(mat_file, check.names = FALSE, stringsAsFactors = FALSE)
kv("matrix", paste(nrow(m), "genes x", ncol(m), "columns"))
count_cols <- grep("_count$", colnames(m), value = TRUE)
counts <- as.matrix(m[, count_cols]); rownames(counts) <- m$gene_id
colnames(counts) <- sub("_count$", "", count_cols)
counts <- round(counts); storage.mode(counts) <- "integer"
annot <- data.frame(ensembl_id = m$gene_id, symbol = m$gene_name, biotype = m$gene_biotype, stringsAsFactors = FALSE)
write.csv(annot, file.path(OUT, "gene_annotation.csv"), row.names = FALSE)

sample_info <- read.csv(file.path(META, "dataset2_samples.csv"), stringsAsFactors = FALSE)
stopifnot(setequal(sample_info$sample_id, colnames(counts)))
counts <- counts[, sample_info$sample_id]
sample_info$genotype  <- factor(sample_info$genotype, levels = c("EK", "pd1245"))
sample_info$diet      <- factor(sample_info$diet, levels = c("HP", "LP"))
sample_info$condition <- factor(paste(sample_info$genotype, sample_info$diet, sep = "_"),
                                levels = c("EK_HP", "EK_LP", "pd1245_HP", "pd1245_LP"))
rownames(sample_info) <- sample_info$sample_id
print(table(sample_info$genotype, sample_info$diet))
check_expected("d2_samples", ncol(counts))
kv("library sizes (M assigned reads)", paste(round(colSums(counts) / 1e6, 1), collapse = " "))
write.csv(counts, file.path(OUT, "raw_counts.csv"), quote = FALSE)
write.csv(sample_info, file.path(OUT, "sample_metadata.csv"), row.names = FALSE)

banner("2. Filter (>=10 counts in >=3 samples)")
keep <- rowSums(counts >= 10) >= 3
kv("genes before / after", c(nrow(counts), sum(keep)))
counts_f <- counts[keep, ]

banner("3. Full model -> VST, PCA (Fig S7D)")
dds <- DESeqDataSetFromMatrix(counts_f, sample_info, design = ~ genotype + diet + genotype:diet)
dds <- DESeq(dds, quiet = TRUE)
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
write.csv(data.frame(gene_id = rownames(vst_mat), gene_name = annot$symbol[match(rownames(vst_mat), annot$ensembl_id)],
                     vst_mat, check.names = FALSE), file.path(OUT, "vst_transformed_counts.csv"), row.names = FALSE)
kv("VST matrix", paste(dim(vst_mat), collapse = " x "))

pcaData <- plotPCA(vsd, intgroup = c("genotype", "diet"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
kv("PC1 / PC2 % variance", percentVar[1:2])
print(pcaData[order(pcaData$PC1), c("genotype", "diet", "PC1", "PC2")])
write.csv(pcaData, file.path(OUT, "PCA_coordinates.csv"), row.names = TRUE)
p <- ggplot(pcaData, aes(PC1, PC2, color = genotype, shape = diet)) + geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + theme_bw() + scale_color_manual(values = GENO_COLORS[c("EK", "pd1245")]) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))
ggsave(file.path(OUT, "PCA_S7D.pdf"), p, width = 6, height = 5)

banner("4. Group model -> pd1245 vs EK within each diet")
dds_cond <- DESeqDataSetFromMatrix(counts_f, sample_info, design = ~ condition)
dds_cond <- DESeq(dds_cond, quiet = TRUE)
for (diet in c("HP", "LP")) {
  nm  <- paste0(diet, "_pd1245_vs_EK")
  res <- results(dds_cond, contrast = c("condition", paste0("pd1245_", diet), paste0("EK_", diet)), alpha = ALPHA)
  df  <- as.data.frame(res[order(res$padj), ]); df$gene_id <- rownames(df)
  df$gene_name <- annot$symbol[match(df$gene_id, annot$ensembl_id)]
  df <- df[, c("gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  write.csv(df, file.path(OUT, paste0("DESeq2_results_", nm, ".csv")), row.names = FALSE)
  sig <- df[!is.na(df$padj) & df$padj < ALPHA & abs(df$log2FoldChange) > LFC, ]
  write.csv(sig, file.path(OUT, paste0("DESeq2_significant_", nm, ".csv")), row.names = FALSE)
  cat(sprintf("  %s: %d DEGs (%d up / %d down)\n", nm, nrow(sig), sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))
}

banner("5. Output verification")
verify_files(file.path(OUT, c("raw_counts.csv", "gene_annotation.csv", "sample_metadata.csv", "vst_transformed_counts.csv",
                              "PCA_coordinates.csv", "PCA_S7D.pdf",
                              "DESeq2_results_HP_pd1245_vs_EK.csv", "DESeq2_results_LP_pd1245_vs_EK.csv",
                              "DESeq2_significant_HP_pd1245_vs_EK.csv", "DESeq2_significant_LP_pd1245_vs_EK.csv")))
write_session()
cat("\nDone: 04_dataset2_DESeq2.R\n")
