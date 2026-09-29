# ============================================================================
# 05_figure2.R  --  Figure 2 (publication panels) from the results of 01 and 03
#   2A : sample x gene z-score heatmap of the 2,279 DEGs (union of the nine
#        contrasts in 2C); DESeq2 normalised counts, per-gene z-score across the
#        18 samples, genes clustered by correlation distance (1 - r)/2, complete
#        linkage; samples ordered by genotype x diet (not clustered).
#   2B : PCA, plotPCA method (top 500 variable VST genes)
#   2C : DEG summary bar chart (DEG_summary_table.csv)
#   2D : WGCNA module-trait correlations with significance stars (grey dropped)
#   2E : five significant KEGG pathways per trait-significant module, chosen by
#        gene Count; dot size = GeneRatio, colour = BH-adjusted p
# Layout 7.5 x 8.4 in, Arial 8-12 pt; single panels + assembled figure to
# figures/Figure2/ (PDF + TIFF).
# ============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), "..", "R", "common.R"))
source(file.path(ROOT, "R", "fig_style.R"))
need_pkgs("ComplexHeatmap", "circlize", "dplyr", "matrixStats", "gplots", "WGCNA")
OUT <- file.path(FIGURES, "Figure2"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ht_opt(message = FALSE, HEATMAP_LEGEND_PADDING = unit(1, "mm"))

LETTER_MARGIN <- 0.18
FIG_W <- 7.5
W_TOP <- (FIG_W - 2 * LETTER_MARGIN) / 2
W_BD  <- round(W_TOP * 2 / 3, 2)
W_E   <- FIG_W - 2 * LETTER_MARGIN - W_BD
H <- c(A = 2.9, C = 2.9, B = 2.4, D = 3.1); H[["E"]] <- H[["B"]] + H[["D"]]
W <- c(A = W_TOP, C = W_TOP, B = W_BD, D = 2.8, E = W_E)
D_X <- 0.1

# ============================================================================
# 1. DATA
# ============================================================================
banner("1. Data")
sample_info <- read.csv(file.path(RES_D1, "sample_metadata.csv"))
stopifnot(setequal(sample_info$condition, COND_LEVELS))
deg_union <- readLines(file.path(RES_D1, "DEG_union_9_comparisons.txt"))
kv("DEG union (2A)", length(deg_union))
norm_counts <- as.matrix(read.csv(file.path(RES_D1, "normalized_counts.csv"), row.names = 1, check.names = FALSE))
stopifnot(all(deg_union %in% rownames(norm_counts)))
mat_A <- norm_counts[deg_union, sample_info$sample_id]
kv("2A matrix", paste(dim(mat_A), collapse = " x "))

# ============================================================================
# 2. PANEL 2A
# ============================================================================
make_panel_A <- function(mat) {
  z <- t(scale(t(mat)))
  hc_genes <- hclust(as.dist((1 - cor(t(mat))) / 2), method = "complete")
  zmax <- 4                                  # max |z| with n = 18 is (n-1)/sqrt(n) = 4.0
  col_fun <- colorRamp2(seq(-zmax, zmax, length.out = 75), gplots::bluered(75))
  cond <- factor(sample_info$condition, levels = COND_LEVELS)
  ht <- Heatmap(t(z), name = "zscore", col = col_fun, cluster_rows = FALSE, row_split = cond, row_title = NULL,
                row_gap = unit(0.6, "mm"), cluster_columns = hc_genes, column_dend_height = unit(0.55, "in"),
                column_dend_gp = gpar(lwd = 0.25), show_row_names = FALSE, show_column_names = FALSE,
                show_heatmap_legend = FALSE,
                left_annotation = rowAnnotation(labels = anno_empty(border = FALSE, width = unit(0.85, "in"))),
                use_raster = TRUE, raster_quality = 4, raster_device = "agg_png")
  lgd <- Legend(col_fun = col_fun, title = "Z-score", direction = "horizontal", at = c(-4, -2, 0, 2, 4),
                legend_width = unit(0.8, "in"), grid_height = unit(2.2, "mm"), title_position = "topcenter",
                title_gp = gpar(fontsize = FS_SMALL, fontface = "bold", fontfamily = FONT),
                labels_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT))
  function() {
    vp <- paste0("panelA_", as.integer(runif(1, 1, 1e9)))
    pushViewport(viewport(name = vp, gp = gpar(fontfamily = FONT)))
    draw(ht, newpage = FALSE, padding = unit(c(1, 1, 1, 1), "mm"))
    for (i in seq_along(COND_LEVELS)) {
      geno <- sub("_.*", "", COND_LEVELS[i]); diet <- sub(".*_", "", COND_LEVELS[i])
      lab <- if (geno == "EK") paste0("EK_", diet) else parse(text = sprintf('paste(italic("%s"), "_%s")', geno, diet))[[1]]
      decorate_annotation("labels", slice = i, {
        grid.lines(x = unit(c(1, 1), "npc") - unit(1, "mm"), y = unit(c(0.06, 0.94), "npc"),
                   gp = gpar(col = GENO_COLORS[geno], lwd = 2.5, lineend = "butt"))
        grid.text(lab, x = unit(1, "npc") - unit(2.5, "mm"), just = "right",
                  gp = gpar(col = GENO_TEXT[geno], fontsize = FS_MAIN, fontfamily = FONT))
      })
    }
    seekViewport(vp)
    draw(lgd, x = unit(1.5, "mm"), y = unit(1, "npc") - unit(1.5, "mm"), just = c("left", "top"))
    popViewport()
  }
}

# ============================================================================
# 3. PANEL 2B
# ============================================================================
vsd_mat <- as.matrix(read.csv(file.path(RES_D1, "vst_transformed_counts.csv"), row.names = 1, check.names = FALSE))[, sample_info$sample_id]
top500 <- order(matrixStats::rowVars(vsd_mat), decreasing = TRUE)[1:500]
pca <- prcomp(t(vsd_mat[top500, ]))
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
pcaData <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], genotype = sample_info$genotype, diet = sample_info$diet)
kv("2B PC1 / PC2 %", percentVar[1:2])
p_B <- ggplot(pcaData, aes(PC1, PC2, color = genotype, shape = diet)) + geom_point(size = 2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + scale_color_manual(values = GENO_COLORS, labels = geno_labels, name = "Genotype") +
  scale_shape_manual(values = c(HP = 16, LP = 17), name = "Diet") + theme_pub() +
  theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"),
        legend.position = "bottom", legend.justification = "center", legend.box = "vertical",
        legend.box.just = "center", legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(2, 4, 2, 4),
        legend.box.background = element_rect(color = "black", linewidth = 0.4, fill = NA),
        legend.box.spacing = unit(2, "mm"), legend.key.size = unit(3, "mm"), panel.grid.minor = element_blank())

# ============================================================================
# 4. PANEL 2C
# ============================================================================
summary_table <- read.csv(file.path(RES_D1, "DEG_summary_table.csv"))
plot_data <- summary_table %>% mutate(Category = ifelse(grepl("HP vs LP", Comparison), "Diet Effects", "Genotype Effects"),
                                      Down_Plot = -Downregulated)
plot_data$Comparison <- factor(plot_data$Comparison, levels = rev(plot_data$Comparison))
p_C <- ggplot(plot_data) +
  geom_col(aes(x = Comparison, y = Down_Plot, fill = "Downregulated"), width = 0.7) +
  geom_col(aes(x = Comparison, y = Upregulated, fill = "Upregulated"), width = 0.7) +
  geom_text(aes(x = Comparison, y = Down_Plot, label = Downregulated), hjust = 1.15, color = "cornflowerblue",
            fontface = "bold", size = FS_SMALL / .pt, family = FONT) +
  geom_text(aes(x = Comparison, y = Upregulated, label = Upregulated), hjust = -0.15, color = "firebrick3",
            fontface = "bold", size = FS_SMALL / .pt, family = FONT) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = function(x) as.expression(lapply(x, italic_pd, bold = TRUE))) +
  scale_y_continuous(labels = abs, limits = c(-950, 580), breaks = c(-800, -400, 0, 400)) +
  scale_fill_manual(values = c("Downregulated" = "cornflowerblue", "Upregulated" = "firebrick3")) +
  coord_flip(clip = "off") +
  labs(title = "Differential Expression Summary", subtitle = "Count of Significant Genes (FDR < 0.05, |LFC| > 1)",
       x = NULL, y = "Number of DEGs") +
  theme_minimal(base_size = FS_SMALL, base_family = FONT) +
  theme(text = element_text(family = FONT, color = "black"), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.text.y = element_text(color = "black", size = FS_SMALL),
        axis.text.x = element_text(color = "black", size = FS_SMALL), axis.title.x = element_text(size = FS_MAIN),
        strip.text = element_text(face = "bold", size = FS_SMALL), strip.background = element_rect(fill = "grey90", color = NA),
        plot.title = element_text(face = "bold", size = FS_MAIN, hjust = 0.5),
        plot.subtitle = element_text(size = FS_SMALL, hjust = 0.5, margin = margin(b = 4)),
        plot.title.position = "plot", legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = FS_SMALL, face = "bold"), legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(3, 3, 3, 3))

# ============================================================================
# 5. PANEL 2D
# ============================================================================
banner("2D module-trait")
mt_cor <- as.matrix(read.csv(file.path(RES_D1_WG, "module_trait_correlations.csv"), row.names = 1))
mt_p   <- as.matrix(read.csv(file.path(RES_D1_WG, "module_trait_pvalues.csv"), row.names = 1))
stopifnot(identical(dim(mt_cor), dim(mt_p)), sum(rownames(mt_cor) == "MEgrey") == 1)
mt_cor <- mt_cor[rownames(mt_cor) != "MEgrey", , drop = FALSE]; mt_p <- mt_p[rownames(mt_p) != "MEgrey", , drop = FALSE]
stars <- ifelse(mt_p < 0.001, "***", ifelse(mt_p < 0.01, "**", ifelse(mt_p < 0.05, "*", "")))
kv("modules shown", nrow(mt_cor)); kv("starred cells", sum(stars != ""))
mod_cols <- setNames(sub("^ME", "", rownames(mt_cor)), rownames(mt_cor))
ht_D <- Heatmap(mt_cor, name = "mt", col = colorRamp2(seq(-1, 1, length.out = 50), WGCNA::blueWhiteRed(50)),
                cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = NA), border = "black",
                row_names_side = "left", row_names_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT),
                column_names_rot = 90, column_names_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT),
                left_annotation = rowAnnotation(module = rownames(mt_cor), col = list(module = mod_cols), show_legend = FALSE,
                                                show_annotation_name = FALSE, simple_anno_size = unit(2.5, "mm")),
                cell_fun = function(j, i, x, y, w, h, fill) if (stars[i, j] != "")
                  grid.text(stars[i, j], x, y - unit(0.6, "mm"), gp = gpar(fontsize = FS_MAIN, fontface = "bold", fontfamily = FONT)),
                heatmap_legend_param = list(title = NULL, at = c(-1, -0.5, 0, 0.5, 1), legend_height = unit(1.3, "in"),
                                            grid_width = unit(2.5, "mm"), labels_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT)))
draw_D <- function() {
  pushViewport(viewport(gp = gpar(fontfamily = FONT)))
  draw(ht_D, newpage = FALSE, padding = unit(c(1, 0, 0.5 * 25.4, 4), "mm"))
  grid.text("Significant Module-Trait\nRelationships", y = unit(1, "npc") - unit(1, "mm"), just = "top",
            gp = gpar(fontsize = FS_MAIN, fontface = "bold", lineheight = 0.95))
  grid.text("(* p<0.05, ** p<0.01, *** p<0.001)", y = unit(1, "npc") - unit(0.36, "in"), just = "top",
            gp = gpar(fontsize = FS_SMALL, fontface = "bold"))
  popViewport()
}

# ============================================================================
# 6. PANEL 2E
# ============================================================================
banner("2E KEGG per module")
# PLLP_USE_PUBLISHED=1 draws the 35 dots of the published panel (Feb 2026 KEGG
# release) from published_results/ instead of this run's KEGG results.
USE_PUBLISHED <- nzchar(Sys.getenv("PLLP_USE_PUBLISHED"))
src_2E <- if (USE_PUBLISHED) file.path(ROOT, "published_results") else RES_D1_WG
kv("2E source", src_2E)
kegg_top5 <- read.csv(file.path(src_2E, "Fig2E_KEGG_top5_per_module.csv"), stringsAsFactors = FALSE)
kegg_top5$Cluster <- tolower(kegg_top5$Cluster)
pref_order <- c("yellow", "greenyellow", "red", "green", "black", "turquoise", "blue")      # published order
mods <- c(pref_order[pref_order %in% kegg_top5$Cluster], setdiff(unique(kegg_top5$Cluster), pref_order))
pretty <- function(m) ifelse(m == "greenyellow", "GreenYellow", paste0(toupper(substr(m, 1, 1)), substr(m, 2, nchar(m))))
kegg_top5$Cluster <- factor(pretty(kegg_top5$Cluster), levels = pretty(mods))
kegg_top5 <- kegg_top5[order(kegg_top5$Cluster), ]
kegg_top5$Description <- factor(kegg_top5$Description, levels = rev(unique(kegg_top5$Description)))
kv("2E dots / pathways / modules", c(nrow(kegg_top5), nlevels(kegg_top5$Description), nlevels(kegg_top5$Cluster)))
write.csv(kegg_top5, file.path(OUT, "Fig2E_KEGG_top5_per_module.csv"), row.names = FALSE)
p_E <- ggplot(kegg_top5, aes(x = Cluster, y = Description, size = GeneRatioNum, color = p.adjust)) + geom_point() +
  scale_color_gradient(low = "red", high = "blue", name = "p.adjust",
                       guide = guide_colorbar(order = 2, barwidth = unit(2.5, "mm"), barheight = unit(0.9, "in"))) +
  scale_size(range = c(0.8, 3), name = "GeneRatio", breaks = c(0.1, 0.2), guide = guide_legend(order = 1)) +
  labs(title = "Five Enriched Pathways per Module", x = NULL, y = NULL) + theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"), plot.title.position = "plot",
        legend.position = "right", legend.key.size = unit(3, "mm"), legend.box.spacing = unit(2, "mm"),
        legend.margin = margin(0, 0, 0, 0))

# ============================================================================
# 7. SAVE
# ============================================================================
banner("Save")
draw_A <- make_panel_A(mat_A); draw_B <- as_draw(p_B); draw_C <- as_draw(p_C); draw_E <- as_draw(p_E)
save_both(draw_A, "Fig2A_allDEGs",      W[["A"]], H[["A"]], OUT)
save_both(draw_B, "Fig2B_PCA",          W[["B"]], H[["B"]], OUT)
save_both(draw_C, "Fig2C_DEG_summary",  W[["C"]], H[["C"]], OUT)
save_both(draw_D, "Fig2D_module_trait", W[["D"]], H[["D"]], OUT)
save_both(draw_E, "Fig2E_KEGG_modules", W[["E"]], H[["E"]], OUT)

FIG_H <- H[["A"]] + H[["E"]]
assemble <- function() {
  panels <- list(list("A", draw_A, 0, 0, NULL), list("C", draw_C, LETTER_MARGIN + W_TOP, 0, NULL),
                 list("B", draw_B, 0, H[["A"]], NULL), list("D", draw_D, 0, H[["A"]] + H[["B"]], D_X),
                 list("E", draw_E, LETTER_MARGIN + W_BD, H[["A"]], NULL))
  for (p in panels) {
    cx <- if (is.null(p[[5]])) p[[3]] + LETTER_MARGIN else p[[5]]
    pushViewport(viewport(x = unit(cx, "in"), y = unit(FIG_H - p[[4]], "in"), width = unit(W[[p[[1]]]], "in"),
                          height = unit(H[[p[[1]]]], "in"), just = c("left", "top")))
    p[[2]](); popViewport()
    grid.text(p[[1]], x = unit(p[[3]] + 0.02, "in"), y = unit(FIG_H - p[[4]] - 0.02, "in"), just = c("left", "top"),
              gp = gpar(fontsize = FS_LETTER, fontface = "bold", fontfamily = FONT))
  }
}
save_both(assemble, "Figure2_full", FIG_W, FIG_H, OUT)
verify_files(file.path(OUT, c(outer(c("Fig2A_allDEGs", "Fig2B_PCA", "Fig2C_DEG_summary", "Fig2D_module_trait",
                                       "Fig2E_KEGG_modules", "Figure2_full"), c(".pdf", ".tiff"), paste0))))
write_session()
cat("\nDone: 05_figure2.R\n")
