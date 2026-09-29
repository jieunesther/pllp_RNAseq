# ============================================================================
# 10_figureS7.R  --  Figure S7D/E (dataset 2, from 04; dataset 1 from 01)
#   S7D : PCA of dataset 2 (EK vs adapted F8 pd1245, HP / LP), plotPCA method
#   S7E : immune-gene heatmaps. Left: dataset 1 EK vs F4 pd1245 (pd1116 samples
#         excluded; z-score over the 12 EK + pd1245 samples). Right: dataset 2
#         EK vs F8 pd1245 (z-score over 14 samples). Group means per genotype x
#         diet. Only curated immune genes that pass the expression filter in
#         dataset 2 are shown (the others have unreliably low counts in the
#         3'-biased library). Markers: * padj < 0.05, ** < 0.01, *** < 0.001 for
#         pd1245 vs EK within diet, on both cells of the pair.
# ============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), "..", "R", "common.R"))
source(file.path(ROOT, "R", "fig_style.R"))
need_pkgs("ComplexHeatmap", "circlize")
ht_opt(message = FALSE)
OUT <- file.path(FIGURES, "FigureS7"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
col_ord <- c("EK_HP", "pd1245_HP", "EK_LP", "pd1245_LP")
tiers <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

# ---- S7D ------------------------------------------------------------------
banner("S7D PCA")
pcaData <- read.csv(file.path(RES_D2, "PCA_coordinates.csv"), row.names = 1)
pv <- round(100 * c(var(pcaData$PC1), var(pcaData$PC2)) / sum(apply(read.csv(file.path(RES_D2, "PCA_coordinates.csv"), row.names = 1)[, c("PC1", "PC2")], 2, var)))
# percentVar as computed by plotPCA is not stored; recompute from the VST matrix
vst2 <- read.csv(file.path(RES_D2, "vst_transformed_counts.csv"), check.names = FALSE)
m2 <- as.matrix(vst2[, -(1:2)]); rownames(m2) <- vst2$gene_id
top500 <- order(apply(m2, 1, var), decreasing = TRUE)[1:500]
pca <- prcomp(t(m2[top500, ])); percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
kv("PC1 / PC2 %", percentVar[1:2]); check_expected("d2_PCA_PC1_pct", percentVar[1]); check_expected("d2_PCA_PC2_pct", percentVar[2])
pd <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], genotype = pcaData[rownames(pca$x), "genotype"], diet = pcaData[rownames(pca$x), "diet"])
p_D <- ggplot(pd, aes(PC1, PC2, color = genotype, shape = diet)) + geom_point(size = 2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() +
  scale_color_manual(values = GENO_COLORS[c("EK", "pd1245")], labels = geno_labels, name = "Genotype") +
  scale_shape_manual(values = c(HP = 16, LP = 17), name = "Diet") + theme_pub() +
  theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.position = "bottom",
        legend.box = "vertical", legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0), legend.key.size = unit(3, "mm"),
        panel.grid.minor = element_blank())
save_both(as_draw(p_D), "FigS7D_PCA", 3.0, 3.0, OUT)

# ---- S7E: dataset 2 panel ---------------------------------------------------
banner("S7E dataset 2 (EK vs F8 pd1245)")
cur <- read.delim(file.path(CURATED, "fig5G_immune.tsv"), stringsAsFactors = FALSE); cur$group <- factor(cur$group, levels = unique(cur$group))
s2 <- read.csv(file.path(RES_D2, "sample_metadata.csv")); stopifnot(identical(s2$sample_id, colnames(m2)))
grp2 <- factor(s2$condition, levels = col_ord)
cur$in_analysis <- cur$ensembl_gene_id %in% rownames(m2)
kv("curated genes detected in dataset 2", sum(cur$in_analysis)); kv("dropped", paste(cur$query_symbol[!cur$in_analysis], collapse = ", "))
check_expected("figS7E_genes", sum(cur$in_analysis))
curP <- cur[cur$in_analysis, ]; curP$group <- droplevels(curP$group)
z2 <- t(scale(t(m2[curP$ensembl_gene_id, ]))); rownames(z2) <- curP$query_symbol
z2g <- sapply(col_ord, function(g) rowMeans(z2[, grp2 == g, drop = FALSE]))
padj2 <- function(d) { x <- read.csv(file.path(RES_D2, paste0("DESeq2_results_", d, "_pd1245_vs_EK.csv"))); setNames(x$padj, x$gene_id)[curP$ensembl_gene_id] }
hp2 <- tiers(padj2("HP")); lp2 <- tiers(padj2("LP"))
disp2 <- cbind(EK_HP = hp2, pd1245_HP = hp2, EK_LP = lp2, pd1245_LP = lp2); rownames(disp2) <- curP$query_symbol
kv("dataset 2 significant (HP / LP)", c(sum(hp2 != ""), sum(lp2 != "")))
check_expected("d2_curated_immune_sig", sum(hp2 != "") + sum(lp2 != ""))

# ---- S7E: dataset 1 panel (EK + pd1245 only) ---------------------------------
banner("S7E dataset 1 (EK vs F4 pd1245, pd1116 excluded)")
vst1 <- as.matrix(read.csv(file.path(RES_D1, "vst_transformed_counts.csv"), row.names = 1, check.names = FALSE))
annot1 <- read.csv(file.path(RES_D1, "gene_annotation.csv"), stringsAsFactors = FALSE)
s1 <- read.csv(file.path(RES_D1, "sample_metadata.csv")); keep <- s1$genotype %in% c("EK", "pd1245")
grp1 <- factor(s1$condition[keep], levels = col_ord)
row1 <- annot1$symbol_used[match(curP$ensembl_gene_id, annot1$ensembl_id)]; stopifnot(!any(is.na(row1)))
z1 <- t(scale(t(vst1[row1, keep]))); rownames(z1) <- curP$query_symbol
z1g <- sapply(col_ord, function(g) rowMeans(z1[, grp1 == g, drop = FALSE]))
padj1 <- function(d) { x <- read.csv(file.path(RES_D1, paste0("DESeq2_results_", d, "_1245_vs_EK.csv"))); setNames(x$padj, x$gene_id)[row1] }
hp1 <- tiers(padj1("HP")); lp1 <- tiers(padj1("LP"))
disp1 <- cbind(EK_HP = hp1, pd1245_HP = hp1, EK_LP = lp1, pd1245_LP = lp1); rownames(disp1) <- curP$query_symbol
kv("dataset 1 significant (HP / LP)", c(sum(hp1 != ""), sum(lp1 != "")))
meta <- data.frame(group = curP$group, gene = curP$query_symbol, ensembl_gene_id = curP$ensembl_gene_id)
write.csv(data.frame(meta, z1g, check.names = FALSE), file.path(RES_HM, "figS7E_dataset1_EKpd1245_zscore_group_matrix.csv"), row.names = FALSE)
write.csv(data.frame(meta, z2g, check.names = FALSE), file.path(RES_HM, "figS7E_dataset2_zscore_group_matrix.csv"), row.names = FALSE)
write.csv(data.frame(meta, d1_HP_padj = padj1("HP"), d1_LP_padj = padj1("LP"), d2_HP_padj = padj2("HP"), d2_LP_padj = padj2("LP")),
          file.path(RES_HM, "figS7E_significance.csv"), row.names = FALSE)

# ---- draw side by side --------------------------------------------------------
lim <- max(abs(c(z1g, z2g))); col_fun <- colorRamp2(c(-lim, 0, lim), c("#2166AC", "#F7F7F7", "#B2182B"))
brk <- pretty(c(-lim, lim), 6); brk <- brk[abs(brk) <= lim]
cat_colors <- c("pro-inflammatory" = "#FF9289", "antimicrobial" = "#FF8AFF", "immune-homeostasis/anti-inflammatory" = "#00DB98",
                "antiviral/defense-response" = "#00CBFF", "T-cell-associated" = "#BEC100")
mk <- function(z, disp, title, show_names, legend) Heatmap(
  z, name = "Z-score", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, column_title = title,
  column_title_gp = gpar(fontsize = FS_MAIN, fontface = "bold", fontfamily = FONT),
  row_split = curP$group, row_title = NULL, row_gap = unit(0.7, "mm"),
  column_split = factor(sub(".*_", "", colnames(z)), levels = c("HP", "LP")), column_gap = unit(1.2, "mm"),
  show_column_names = FALSE, show_row_names = show_names, row_names_side = "right",
  row_names_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT), rect_gp = gpar(col = "grey70", lwd = 0.3),
  top_annotation = HeatmapAnnotation(genotype = sub("_(HP|LP)$", "", colnames(z)), col = list(genotype = GENO_COLORS),
                                     simple_anno_size = unit(2.2, "mm"), show_legend = FALSE, show_annotation_name = FALSE),
  left_annotation = if (!show_names) rowAnnotation(category = curP$group, col = list(category = cat_colors), simple_anno_size = unit(2.2, "mm"),
                                                   show_legend = FALSE, show_annotation_name = FALSE) else NULL,
  cell_fun = function(j, i, x, y, w, h, fill) if (disp[i, j] != "") grid.text(disp[i, j], x, y, gp = gpar(fontsize = FS_MARK, fontfamily = FONT)),
  show_heatmap_legend = legend,
  heatmap_legend_param = list(title = "Z-score", at = brk, title_gp = gpar(fontsize = FS_SMALL, fontface = "bold", fontfamily = FONT),
                              labels_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT), grid_width = unit(2.5, "mm"), legend_height = unit(1.1, "in")))
ht1 <- mk(z1g, disp1, expression(paste("EK vs F4 ", italic("pd1245"))), FALSE, FALSE)
ht2 <- mk(z2g, disp2, expression(paste("EK vs F8 ", italic("pd1245"))), TRUE, TRUE)
lgd_cat <- Legend(labels = c("pro-inflammatory", "antimicrobial", "immune-homeostasis/\nanti-inflammatory", "antiviral/\ndefense-response", "T-cell-associated")[
  match(levels(curP$group), names(cat_colors))], title = "gene category", legend_gp = gpar(fill = unname(cat_colors[levels(curP$group)])),
  title_gp = gpar(fontsize = FS_SMALL, fontface = "bold", fontfamily = FONT), labels_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT),
  grid_width = unit(2.5, "mm"), grid_height = unit(2.5, "mm"), row_gap = unit(0.8, "mm"))
draw_E <- function() {
  pushViewport(viewport(gp = gpar(fontfamily = FONT)))
  pushViewport(viewport(y = unit(1, "npc"), height = unit(0.32, "in"), just = "top"))
  grid.text(expression(paste("*, **, *** : ", italic("p"), "-adj < 0.05, 0.01, 0.001,  ", italic("pd1245"), " vs EK within diet")),
            x = unit(2.5, "mm"), just = "left", gp = gpar(fontsize = FS_SMALL, fontfamily = FONT))
  popViewport()
  pushViewport(viewport(y = 0, height = unit(1, "npc") - unit(0.32, "in"), just = "bottom"))
  draw(ht1 + ht2, annotation_legend_list = list(lgd_cat), heatmap_legend_side = "right", annotation_legend_side = "right",
       merge_legend = TRUE, newpage = FALSE, ht_gap = unit(4, "mm"), padding = unit(c(1, 1, 1, 1), "mm"))
  popViewport(); popViewport()
}
save_both(draw_E, "FigS7E_immune_heatmaps", 5.2, 0.32 + 0.6 + 0.135 * nrow(z1g) + 0.4, OUT)
verify_files(file.path(OUT, c("FigS7D_PCA.pdf", "FigS7D_PCA.tiff", "FigS7E_immune_heatmaps.pdf", "FigS7E_immune_heatmaps.tiff")))
write_session()
cat("\nDone: 10_figureS7.R\n")
