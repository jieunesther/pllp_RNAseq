# ============================================================================
# 08_figureS2.R  --  Figure S2 (from 01 and 03)
#   S2A : Venn, diet effects (HP vs LP) in EK / pd1116 / pd1245
#   S2B : Venn, DEGs vs EK under HP (pd1116, pd1245)
#   S2C : Venn, DEGs vs EK under LP
#   S2D : module eigengene box plots, 12 modules (grey excluded), 3 replicates
#         per genotype x diet group
# ============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), "..", "R", "common.R"))
source(file.path(ROOT, "R", "fig_style.R"))
need_pkgs("VennDiagram")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
OUT <- file.path(FIGURES, "FigureS2"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIG_W <- 7.5; FIG_H <- 8.6; H_VENN <- 2.30; H_D <- FIG_H - H_VENN; W_VENN <- FIG_W / 3

sig_genes <- function(f) { d <- read.csv(file.path(RES_D1, f)); unique(d$gene_id[!is.na(d$padj) & d$padj < ALPHA & abs(d$log2FoldChange) > LFC]) }
S <- list(EK_diet = sig_genes("DESeq2_significant_EK_HP_vs_LP.csv"), p1116_diet = sig_genes("DESeq2_significant_1116_HP_vs_LP.csv"),
          p1245_diet = sig_genes("DESeq2_significant_1245_HP_vs_LP.csv"), HP_1116 = sig_genes("DESeq2_significant_HP_1116_vs_EK.csv"),
          HP_1245 = sig_genes("DESeq2_significant_HP_1245_vs_EK.csv"), LP_1116 = sig_genes("DESeq2_significant_LP_1116_vs_EK.csv"),
          LP_1245 = sig_genes("DESeq2_significant_LP_1245_vs_EK.csv"))
banner("Venn set sizes"); print(sapply(S, length))

venn_grob <- function(sets, fills, cat_pos, cat_dist) {
  venn.diagram(sets, filename = NULL, fill = fills, alpha = 0.5, lwd = 0.6, scaled = FALSE, cex = 1, cat.cex = 1,
               fontface = "bold", cat.fontface = "bold", fontfamily = FONT, cat.fontfamily = FONT, cat.col = "red",
               cat.pos = cat_pos, cat.dist = cat_dist, print.mode = c("raw", "percent"), margin = 0.14)
}
g_A <- venn_grob(setNames(S[c("EK_diet", "p1116_diet", "p1245_diet")], c("EK\nHP vs LP", "pd1116\nHP vs LP", "pd1245\nHP vs LP")),
                 unname(GENO_COLORS), cat_pos = c(-27, 27, 180), cat_dist = c(0.055, 0.055, 0.105))
# Two-set panels are drawn directly so pd1116 (green) is always on the left:
# VennDiagram would place the larger set on the left, swapping B and C.
pairwise_venn <- function(nA, nB, nAB, labA, labB, fillA = GENO_COLORS[["pd1116"]], fillB = GENO_COLORS[["pd1245"]]) function() {
  total <- nA + nB - nAB; lab_n <- function(k) sprintf("%d\n(%.3g%%)", k, 100 * k / total)
  pushViewport(viewport(width = unit(1, "snpc"), height = unit(1, "snpc"), gp = gpar(fontsize = FS_SMALL, fontfamily = FONT)))
  r <- 0.26; cy <- 0.42; cxA <- 0.33; cxB <- 0.67
  grid.circle(cxA, cy, r, gp = gpar(fill = adjustcolor(fillA, alpha.f = 0.5), col = "black", lwd = 0.8))
  grid.circle(cxB, cy, r, gp = gpar(fill = adjustcolor(fillB, alpha.f = 0.5), col = "black", lwd = 0.8))
  grid.text(lab_n(nA - nAB), cxA - 0.52 * r, cy, gp = gpar(fontface = "bold"))
  grid.text(lab_n(nAB), 0.5, cy, gp = gpar(fontface = "bold"))
  grid.text(lab_n(nB - nAB), cxB + 0.52 * r, cy, gp = gpar(fontface = "bold"))
  grid.text(labA, cxA - 0.05, cy + r + 0.05, just = c("centre", "bottom"), gp = gpar(col = "red", fontface = "bold"))
  grid.text(labB, cxB + 0.05, cy + r + 0.05, just = c("centre", "bottom"), gp = gpar(col = "red", fontface = "bold"))
  popViewport()
}
draw_B <- pairwise_venn(length(S$HP_1116), length(S$HP_1245), length(intersect(S$HP_1116, S$HP_1245)), "HP\npd1116 vs EK", "HP\npd1245 vs EK")
draw_C <- pairwise_venn(length(S$LP_1116), length(S$LP_1245), length(intersect(S$LP_1116, S$LP_1245)), "LP\npd1116 vs EK", "LP\npd1245 vs EK")
draw_venn <- function(g) function() { pushViewport(viewport(gp = gpar(fontsize = FS_SMALL, fontfamily = FONT))); grid.draw(g); popViewport() }

# ---- S2D: eigengenes saved by 03 --------------------------------------------
banner("S2D eigengenes")
MEs <- read.csv(file.path(RES_D1_WG, "module_eigengenes.csv"), check.names = FALSE)
sample_info <- read.csv(file.path(RES_D1, "sample_metadata.csv"))
stopifnot(identical(MEs$sample_id, sample_info$sample_id))
mt <- read.csv(file.path(RES_D1_WG, "module_trait_correlations.csv"), row.names = 1)
mod_order <- sub("^ME", "", rownames(mt)); mod_order <- mod_order[mod_order != "grey"]
kv("modules", paste(mod_order, collapse = ", "))
d_long <- do.call(rbind, lapply(mod_order, function(m) data.frame(
  module = factor(paste("Module", m), levels = paste("Module", mod_order)), color = m,
  condition = factor(sample_info$condition, levels = COND_LEVELS), value = MEs[[paste0("ME", m)]])))
p_D <- ggplot(d_long, aes(x = condition, y = value, fill = module)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.25, fatten = 1.4, width = 0.65) +
  facet_wrap(~ module, nrow = 3, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = setNames(mod_order, paste("Module", mod_order)), guide = "none") +
  scale_x_discrete(labels = function(x) as.expression(lapply(x, italic_pd))) +
  labs(x = NULL, y = "Eigengene Expression") + theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = FS_SMALL), axis.text.y = element_text(size = FS_SMALL),
        axis.title.y = element_text(size = FS_MAIN), strip.text = element_text(size = FS_MAIN, face = "bold"),
        strip.background = element_blank(), panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4), panel.spacing = unit(2, "mm"), plot.margin = margin(3, 3, 3, 3))

save_both(draw_venn(g_A), "FigS2A_venn_diet", W_VENN, H_VENN, OUT)
save_both(draw_B, "FigS2B_venn_HP", W_VENN, H_VENN, OUT)
save_both(draw_C, "FigS2C_venn_LP", W_VENN, H_VENN, OUT)
save_both(as_draw(p_D), "FigS2D_module_eigengenes", FIG_W, H_D, OUT)
save_both(function() {
  place_panel("A", draw_venn(g_A), x = 0, top = 0, w = W_VENN, h = H_VENN, fig_h = FIG_H, inset = 0.10)
  place_panel("B", draw_B, x = W_VENN, top = 0, w = W_VENN, h = H_VENN, fig_h = FIG_H, inset = 0.10)
  place_panel("C", draw_C, x = 2 * W_VENN, top = 0, w = W_VENN, h = H_VENN, fig_h = FIG_H, inset = 0.10)
  place_panel("D", as_draw(p_D), x = 0, top = H_VENN, w = FIG_W, h = H_D, fig_h = FIG_H, inset = 0)
}, "FigureS2_full", FIG_W, FIG_H, OUT)
verify_files(file.path(OUT, c(outer(c("FigS2A_venn_diet", "FigS2B_venn_HP", "FigS2C_venn_LP", "FigS2D_module_eigengenes", "FigureS2_full"), c(".pdf", ".tiff"), paste0))))
write_session()
cat("\nDone: 08_figureS2.R\n")
