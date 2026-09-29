# ============================================================================
# 06_figure5AB.R  --  Figure 5A/5B: GO BP enrichment dot plots (from 02)
#   5A : High-protein diet  (Enrichment_GO_HP.csv)
#   5B : Low-protein diet   (Enrichment_GO_LP.csv)
# Columns: pd1116 up / down, pd1245 up / down (vs EK). Rows: GO BP terms that are
# among the 10 most significant (by BH p) in at least one column. Dot size =
# gene count, colour = FDR. Published panels: 11 terms (A), 19 terms (B).
# ============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), "..", "R", "common.R"))
source(file.path(ROOT, "R", "fig_style.R"))
need_pkgs("dplyr")
OUT <- file.path(FIGURES, "Figure5"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

FIG_W <- 7.5; FIG_H <- 3.5; LETTER_MARGIN <- 0.18
W_PANEL <- (FIG_W - 2 * LETTER_MARGIN) / 2; H_PANEL <- FIG_H
LEVELS <- list(HP = c("1116_vs_EK_HP_Up", "1116_vs_EK_HP_Down", "1245_vs_EK_HP_Up", "1245_vs_EK_HP_Down"),
               LP = c("1116_vs_EK_LP_Up", "1116_vs_EK_LP_Down", "1245_vs_EK_LP_Up", "1245_vs_EK_LP_Down"))
# display-only shortening so every label fits on one line at 8 pt
ABBREV <- c("biological process involved in interspecies interaction between organisms" = "interspecies interaction btw. organisms",
            "regulation of receptor signaling pathway via JAK-STAT" = "reg. of receptor signaling via JAK-STAT",
            "regulation of receptor signaling pathway via STAT" = "reg. of receptor signaling via STAT",
            "external encapsulating structure organization" = "external encapsulating structure org.",
            "negative regulation of peptidase activity" = "neg. regulation of peptidase activity",
            "negative regulation of hydrolase activity" = "neg. regulation of hydrolase activity",
            "extracellular structure organization" = "extracellular structure org.",
            "extracellular matrix organization" = "extracellular matrix org.")
short_label <- function(x) ifelse(x %in% names(ABBREV), ABBREV[x], x)

# PLLP_USE_PUBLISHED=1 draws from published_results/ (tables of the Feb 2026 run)
USE_PUBLISHED <- nzchar(Sys.getenv("PLLP_USE_PUBLISHED"))
SHORT2LONG <- c("1116_HP_Up" = "1116_vs_EK_HP_Up", "1116_HP_Down" = "1116_vs_EK_HP_Down",
                "1245_HP_Up" = "1245_vs_EK_HP_Up", "1245_HP_Down" = "1245_vs_EK_HP_Down",
                "1116_LP_Up" = "1116_vs_EK_LP_Up", "1116_LP_Down" = "1116_vs_EK_LP_Down",
                "1245_LP_Up" = "1245_vs_EK_LP_Up", "1245_LP_Down" = "1245_vs_EK_LP_Down")
load_panel <- function(diet) {
  src <- if (USE_PUBLISHED) file.path(ROOT, "published_results") else RES_D1_ENR
  d <- read.csv(file.path(src, paste0("Enrichment_GO_", diet, ".csv")), stringsAsFactors = FALSE)
  d$Cluster <- ifelse(d$Cluster %in% names(SHORT2LONG), SHORT2LONG[d$Cluster], d$Cluster)
  banner(paste("GO BP", diet, if (USE_PUBLISHED) "(published tables)" else "(this run)"))
  kv("source", src); kv("rows", nrow(d)); kv("clusters", unique(d$Cluster))
  stopifnot(all(d$Cluster %in% LEVELS[[diet]]))
  top <- d %>% group_by(Cluster) %>% slice_min(order_by = p.adjust, n = 10) %>% pull(Description)
  pd  <- d %>% filter(Description %in% top)
  kv("terms drawn", length(unique(pd$Description))); print(table(factor(pd$Cluster, levels = LEVELS[[diet]])))
  pd$Cluster <- factor(pd$Cluster, levels = LEVELS[[diet]])
  pd$Description <- factor(pd$Description, levels = sort(unique(pd$Description)))
  pd$Diet <- diet; pd
}
dat_HP <- load_panel("HP"); dat_LP <- load_panel("LP")
check_expected("fig5A_terms", length(unique(dat_HP$Description)))
check_expected("fig5B_terms", length(unique(dat_LP$Description)))
write.csv(rbind(dat_HP, dat_LP), file.path(OUT, "Fig5AB_plotted_terms.csv"), row.names = FALSE)

x_labels <- function(x) as.expression(lapply(x, function(s) {
  geno <- sub("^([0-9]+).*", "pd\\1", s); dir <- if (grepl("_Up$", s)) "up" else "down"
  bquote(paste(italic(.(geno)), " ", .(dir))) }))
make_dotplot <- function(d, title) {
  ggplot(d, aes(x = Cluster, y = Description, size = Count, color = p.adjust)) + geom_point() +
    scale_x_discrete(drop = FALSE, labels = x_labels) + scale_y_discrete(labels = short_label) +
    scale_color_gradient(low = "red", high = "blue", name = "FDR",
                         guide = guide_colorbar(order = 2, barwidth = unit(2.2, "mm"), barheight = unit(0.75, "in"))) +
    scale_size(name = "Count", range = c(0.8, 3.2), guide = guide_legend(order = 1)) +
    labs(title = title, x = NULL, y = NULL) + theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = FS_SMALL, color = "black"),
          axis.text.y = element_text(size = FS_SMALL, color = "black"), plot.title.position = "plot",
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25), panel.grid.minor = element_blank(),
          legend.position = "right", legend.key.size = unit(3, "mm"), legend.box.spacing = unit(1.5, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(1, "mm"), plot.margin = margin(3, 3, 3, 3))
}
p_A <- make_dotplot(dat_HP, "GO BP: High Protein"); p_B <- make_dotplot(dat_LP, "GO BP: Low Protein")
save_both(as_draw(p_A), "Fig5A_GO_HP", W_PANEL, H_PANEL, OUT)
save_both(as_draw(p_B), "Fig5B_GO_LP", W_PANEL, H_PANEL, OUT)
save_both(function() {
  place_panel("A", as_draw(p_A), x = 0, top = 0, w = W_PANEL, h = H_PANEL, fig_h = FIG_H)
  place_panel("B", as_draw(p_B), x = LETTER_MARGIN + W_PANEL, top = 0, w = W_PANEL, h = H_PANEL, fig_h = FIG_H)
}, "Figure5AB", FIG_W, FIG_H, OUT)
verify_files(file.path(OUT, c(outer(c("Fig5A_GO_HP", "Fig5B_GO_LP", "Figure5AB"), c(".pdf", ".tiff"), paste0), "Fig5AB_plotted_terms.csv")))
write_session()
cat("\nDone: 06_figure5AB.R\n")
