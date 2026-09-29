# ============================================================================
# 07_heatmaps_curated_genes.R  --  Fig 3A, 3B, 5G, S4A (dataset 1)
# Group-mean z-score heatmaps of curated gene sets:
#   value   : per-gene z-score of VST across all 18 samples, averaged within
#             each genotype x diet group (6 columns)
#   markers : * = padj < 0.05 EK vs pd1116 (drawn on the EK and pd1116 cells)
#             # = padj < 0.05 pd1116 vs pd1245 (drawn on the pd1116 and pd1245 cells)
#             within the same diet, from the DESeq2 contrasts of script 01
# Genes are resolved by Ensembl ID through results/dataset1/gene_annotation.csv
# (symbol used as a fallback when a curated entry has no Ensembl ID), so the
# panels do not depend on the org.Dr.eg.db symbol version.
# Matrices behind each panel -> results/heatmap_matrices/<panel>_*.csv
# Panels -> figures/Figure3/, figures/Figure5/, figures/FigureS4/
# ============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), "..", "R", "common.R"))
source(file.path(ROOT, "R", "fig_style.R"))
need_pkgs("ComplexHeatmap", "circlize")
ht_opt(message = FALSE)

# ---- dataset-1 inputs -----------------------------------------------------
vst <- as.matrix(read.csv(file.path(RES_D1, "vst_transformed_counts.csv"), row.names = 1, check.names = FALSE))
annot <- read.csv(file.path(RES_D1, "gene_annotation.csv"), stringsAsFactors = FALSE)
sample_info <- read.csv(file.path(RES_D1, "sample_metadata.csv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(vst), sample_info$sample_id))
grp_levels <- c("EK_HP", "pd1116_HP", "pd1245_HP", "EK_LP", "pd1116_LP", "pd1245_LP")
grp <- factor(sample_info$condition, levels = grp_levels)
z_all <- t(scale(t(vst)))                                 # per-gene z across 18 samples
padj_by_row <- function(f) { d <- read.csv(file.path(RES_D1, f)); setNames(d$padj, d$gene_id) }
star_HP <- padj_by_row("DESeq2_results_HP_1116_vs_EK.csv");   star_LP <- padj_by_row("DESeq2_results_LP_1116_vs_EK.csv")
hash_HP <- padj_by_row("DESeq2_results_HP_1245_vs_1116.csv"); hash_LP <- padj_by_row("DESeq2_results_LP_1245_vs_1116.csv")
sig <- function(p) !is.na(p) & p < ALPHA

resolve_rows <- function(cur) {
  by_ens <- annot$symbol_used[match(cur$ensembl_gene_id, annot$ensembl_id)]
  by_sym <- annot$symbol_used[match(tolower(cur$query_symbol), tolower(annot$symbol_used))]
  row <- ifelse(!is.na(by_ens), by_ens, by_sym)
  cur$vst_row <- row
  miss <- is.na(row) | !(row %in% rownames(vst))
  if (any(miss)) cat("  DROPPED (not in filtered VST):", paste(cur$query_symbol[miss], collapse = ", "), "\n")
  cur[!miss, ]
}

build_panel <- function(tag, tsv) {
  banner(tag)
  cur <- read.delim(file.path(CURATED, tsv), stringsAsFactors = FALSE)
  cur$group <- factor(cur$group, levels = unique(cur$group))
  kv("curated genes", nrow(cur)); cur <- resolve_rows(cur); cur$group <- droplevels(cur$group)
  kv("genes plotted", nrow(cur)); print(table(cur$group))
  z_cur <- z_all[cur$vst_row, , drop = FALSE]; rownames(z_cur) <- cur$query_symbol
  z_grp <- sapply(grp_levels, function(g) rowMeans(z_cur[, grp == g, drop = FALSE]))
  id <- cur$vst_row
  sHP <- sig(star_HP[id]); sLP <- sig(star_LP[id]); hHP <- sig(hash_HP[id]); hLP <- sig(hash_LP[id])
  disp <- cbind(EK_HP = ifelse(sHP, "*", ""), pd1116_HP = paste0(ifelse(sHP, "*", ""), ifelse(hHP, "#", "")),
                pd1245_HP = ifelse(hHP, "#", ""), EK_LP = ifelse(sLP, "*", ""),
                pd1116_LP = paste0(ifelse(sLP, "*", ""), ifelse(hLP, "#", "")), pd1245_LP = ifelse(hLP, "#", ""))
  rownames(disp) <- cur$query_symbol
  kv("markers * HP/LP, # HP/LP", c(sum(sHP), sum(sLP), sum(hHP), sum(hLP)))
  kv("z range", sprintf("%.2f", range(z_grp)))
  meta <- data.frame(group = cur$group, gene = cur$query_symbol, ensembl_gene_id = cur$ensembl_gene_id, vst_row = cur$vst_row)
  write.csv(data.frame(meta, z_grp, check.names = FALSE), file.path(RES_HM, paste0(tag, "_zscore_group_matrix.csv")), row.names = FALSE)
  write.csv(data.frame(meta, z_cur, check.names = FALSE), file.path(RES_HM, paste0(tag, "_zscore_persample_matrix.csv")), row.names = FALSE)
  write.csv(data.frame(meta, HP_EKvs1116_padj = star_HP[id], LP_EKvs1116_padj = star_LP[id],
                       HP_1116vs1245_padj = hash_HP[id], LP_1116vs1245_padj = hash_LP[id], disp, check.names = FALSE),
            file.path(RES_HM, paste0(tag, "_significance.csv")), row.names = FALSE)
  list(z = z_grp, disp = disp, groups = cur$group)
}

# ---- generic publication renderer (layout of the published 5G / S4A) -------
draw_group_heatmap <- function(P, cat_colors, cat_title, cat_labels = NULL, key_h = 0.42, geno_h = 0.26) {
  multi <- nlevels(P$groups) > 1
  col_split <- factor(sub(".*_", "", colnames(P$z)), levels = c("HP", "LP"))
  genotype  <- sub("_(HP|LP)$", "", colnames(P$z))
  lim <- max(abs(P$z)); col_fun <- colorRamp2(c(-lim, 0, lim), c("#2166AC", "#F7F7F7", "#B2182B"))
  brk <- pretty(c(-lim, lim), 6); brk <- brk[abs(brk) <= lim]
  ht <- Heatmap(P$z, name = "Z-score", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
                row_split = if (multi) P$groups else NULL, row_title = NULL, row_gap = unit(0.7, "mm"),
                column_split = col_split, column_gap = unit(1.2, "mm"),
                column_title_gp = gpar(fontsize = FS_MAIN, fontface = "bold", fontfamily = FONT),
                show_column_names = FALSE, row_names_side = "right",
                row_names_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT), rect_gp = gpar(col = "grey70", lwd = 0.3),
                top_annotation = HeatmapAnnotation(genotype = genotype, col = list(genotype = GENO_COLORS),
                                                   simple_anno_size = unit(2.2, "mm"), show_legend = FALSE, show_annotation_name = FALSE),
                left_annotation = if (multi) rowAnnotation(category = P$groups, col = list(category = cat_colors),
                                                           simple_anno_size = unit(2.2, "mm"), show_legend = FALSE, show_annotation_name = FALSE) else NULL,
                cell_fun = function(j, i, x, y, w, h, fill) if (P$disp[i, j] != "")
                  grid.text(P$disp[i, j], x, y, gp = gpar(fontsize = FS_MARK, fontfamily = FONT, col = "black")),
                heatmap_legend_param = list(title = "Z-score", at = brk,
                                            title_gp = gpar(fontsize = FS_SMALL, fontface = "bold", fontfamily = FONT),
                                            labels_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT),
                                            grid_width = unit(2.5, "mm"), legend_height = unit(1.1, "in")))
  lgd_cat <- if (multi) Legend(labels = if (is.null(cat_labels)) levels(P$groups) else cat_labels, title = cat_title,
                               legend_gp = gpar(fill = unname(cat_colors[levels(P$groups)])),
                               title_gp = gpar(fontsize = FS_SMALL, fontface = "bold", fontfamily = FONT),
                               labels_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT),
                               grid_width = unit(2.5, "mm"), grid_height = unit(2.5, "mm"), row_gap = unit(0.8, "mm")) else NULL
  function() {
    pushViewport(viewport(gp = gpar(fontfamily = FONT)))
    pushViewport(viewport(y = unit(1, "npc"), height = unit(key_h, "in"), just = "top"))
    grid.rect(width = unit(1, "npc") - unit(2, "mm"), height = unit(1, "npc") - unit(1.5, "mm"), gp = gpar(col = "black", fill = NA, lwd = 0.6))
    grid.text(expression(paste("* : ", italic("p"), "-adj < 0.05, EK vs ", italic("pd1116"))), x = unit(2.5, "mm"), y = unit(0.66, "npc"),
              just = "left", gp = gpar(fontsize = FS_SMALL, fontfamily = FONT))
    grid.text(expression(paste("# : ", italic("p"), "-adj < 0.05, ", italic("pd1116"), " vs ", italic("pd1245"))), x = unit(2.5, "mm"),
              y = unit(0.30, "npc"), just = "left", gp = gpar(fontsize = FS_SMALL, fontfamily = FONT))
    popViewport()
    pushViewport(viewport(y = unit(0, "npc"), height = unit(geno_h, "in"), just = "bottom"))
    x_cur <- unit(2.5, "mm")
    for (g in names(GENO_COLORS)) {
      grid.rect(x = x_cur, y = unit(0.5, "npc"), width = unit(2.2, "mm"), height = unit(2.2, "mm"), just = c("left", "centre"),
                gp = gpar(fill = GENO_COLORS[[g]], col = "grey40", lwd = 0.3))
      lab <- if (g == "EK") "EK" else parse(text = sprintf('italic("%s")', g))[[1]]
      grid.text(lab, x = x_cur + unit(3, "mm"), y = unit(0.5, "npc"), just = c("left", "centre"), gp = gpar(fontsize = FS_SMALL, fontfamily = FONT))
      x_cur <- x_cur + unit(3, "mm") + convertWidth(grobWidth(textGrob(lab, gp = gpar(fontsize = FS_SMALL, fontfamily = FONT))), "in") + unit(4, "mm")
    }
    popViewport()
    pushViewport(viewport(y = unit(geno_h, "in"), just = "bottom", height = unit(1, "npc") - unit(key_h + geno_h, "in")))
    draw(ht, annotation_legend_list = if (multi) list(lgd_cat) else list(), heatmap_legend_side = "right",
         annotation_legend_side = "right", merge_legend = TRUE, newpage = FALSE, padding = unit(c(1, 1, 1, 1), "mm"))
    popViewport(); popViewport()
  }
}

# ---- Fig 5G ----------------------------------------------------------------
P5G <- build_panel("fig5G", "fig5G_immune.tsv")
check_expected("fig5G_genes", nrow(P5G$z))
cat5G <- c("pro-inflammatory" = "#FF9289", "antimicrobial" = "#FF8AFF", "immune-homeostasis/anti-inflammatory" = "#00DB98",
           "antiviral/defense-response" = "#00CBFF", "T-cell-associated" = "#BEC100")
lab5G <- c("pro-inflammatory", "antimicrobial", "immune-homeostasis/\nanti-inflammatory", "antiviral/\ndefense-response", "T-cell-associated")
save_both(draw_group_heatmap(P5G, cat5G, "gene category", lab5G), "Fig5G_immune_heatmap", 3.2, 5.2, file.path(FIGURES, "Figure5"))

# ---- Fig S4A ---------------------------------------------------------------
PS4 <- build_panel("figS4A", "figS4A_celltype_markers.tsv")
check_expected("figS4A_genes", nrow(PS4$z))
catS4 <- setNames(c("#FF9289", "#FF8AFF", "#00DB98", "#00CBFF", "#BEC100"), levels(PS4$groups))
save_both(draw_group_heatmap(PS4, catS4, "cell type"), "FigS4A_celltype_markers", 3.9, 7.6, file.path(FIGURES, "FigureS4"))

# ---- Fig 3A, 3B ------------------------------------------------------------
P3A <- build_panel("fig3A", "fig3A_endocytic_machinery.tsv")
save_both(draw_group_heatmap(P3A, NULL, NULL), "Fig3A_endocytic_heatmap", 3.2, 0.42 + 0.26 + 0.55 + 0.135 * nrow(P3A$z) + 0.3, file.path(FIGURES, "Figure3"))
P3B <- build_panel("fig3B", "fig3B_protease_lysosome.tsv")
save_both(draw_group_heatmap(P3B, NULL, NULL), "Fig3B_protease_lysosome_heatmap", 3.2, 0.42 + 0.26 + 0.55 + 0.135 * nrow(P3B$z) + 0.3, file.path(FIGURES, "Figure3"))

verify_files(c(file.path(FIGURES, "Figure5", c("Fig5G_immune_heatmap.pdf", "Fig5G_immune_heatmap.tiff")),
               file.path(FIGURES, "FigureS4", c("FigS4A_celltype_markers.pdf", "FigS4A_celltype_markers.tiff")),
               file.path(FIGURES, "Figure3", c("Fig3A_endocytic_heatmap.pdf", "Fig3B_protease_lysosome_heatmap.pdf")),
               file.path(RES_HM, paste0(c("fig5G", "figS4A", "fig3A", "fig3B"), "_zscore_group_matrix.csv"))))
write_session()
cat("\nDone: 07_heatmaps_curated_genes.R\n")
