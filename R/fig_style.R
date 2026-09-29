# ============================================================================
# fig_style.R  --  shared publication style for the figure scripts (05-10)
# ============================================================================
# Journal specs enforced here (PLoS): Arial 8-12 pt, panel <= 7.5 x 8.75 in,
# TIFF 600 dpi LZW (falls back to 300 dpi if a file exceeds 10 MB) + PDF.
# On Linux "Arial" resolves to the metric-compatible Liberation Sans.
# ============================================================================
need_pkgs("ggplot2", "grid", "ragg")

FONT      <- "Arial"
FS_SMALL  <- 8    # tick labels, legend text
FS_MAIN   <- 9    # axis / legend titles, panel titles
FS_LETTER <- 12   # panel letters
FS_MARK   <- 6.5  # * / # heatmap markers (graphical symbols, intentionally < 8 pt)

# grid measures text on a pdf() null device that only knows the PostScript
# families, so it warns "font family 'Arial' not found". Rendering goes through
# cairo / ragg, which resolve Arial correctly; muffle only that one message.
muffle_font_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("not found in PostScript font database", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}

theme_pub <- function() {
  theme_bw(base_size = FS_SMALL, base_family = FONT) +
    theme(text          = element_text(family = FONT, size = FS_SMALL, color = "black"),
          axis.text     = element_text(size = FS_SMALL, color = "black"),
          axis.title    = element_text(size = FS_MAIN),
          legend.title  = element_text(size = FS_MAIN),
          legend.text   = element_text(size = FS_SMALL),
          plot.title    = element_text(size = FS_MAIN, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = FS_SMALL, hjust = 0.5),
          plot.margin   = margin(3, 3, 3, 3))
}

# "HP: pd1116 vs EK" -> plotmath with pd alleles in italics
italic_pd <- function(s, bold = FALSE) {
  if (bold) {
    s2 <- gsub("(pd[0-9]+)", '"), bolditalic("\\1"), bold("', s)
    parse(text = paste0('paste(bold("', s2, '"))'))[[1]]
  } else {
    s2 <- gsub("(pd[0-9]+)", '", italic("\\1"), "', s)
    parse(text = paste0('paste("', s2, '")'))[[1]]
  }
}
geno_labels <- function(x) as.expression(lapply(x, function(s) if (grepl("^pd", s)) bquote(italic(.(s))) else s))

open_pdf <- if (capabilities("aqua")) {
  function(file, width, height) grDevices::quartz(type = "pdf", file = file, width = width, height = height, family = FONT)
} else if (capabilities("cairo")) {
  function(file, width, height) grDevices::cairo_pdf(file, width = width, height = height, family = FONT)
} else stop("No quartz or cairo PDF device available.")

# draw_fun: function of no arguments drawing on the current device.
save_both <- function(draw_fun, name, w, h, output_dir) muffle_font_warning({
  stopifnot(w <= 7.5 + 1e-9, h <= 8.75 + 1e-9)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_f <- file.path(output_dir, paste0(name, ".pdf"))
  open_pdf(pdf_f, w, h); grid.newpage(); draw_fun(); invisible(dev.off())
  tif_f <- file.path(output_dir, paste0(name, ".tiff"))
  for (res in c(600, 300)) {
    ragg::agg_tiff(tif_f, width = w, height = h, units = "in", res = res, compression = "lzw", background = "white")
    grid.newpage(); draw_fun(); invisible(dev.off())
    if (file.size(tif_f) < 10e6) break
  }
  cat(sprintf("Saved %-28s %.2f x %.2f in | TIFF %d dpi, %.1f MB\n", name, w, h, res, file.size(tif_f) / 1e6))
})
as_draw <- function(p) function() grid.draw(ggplotGrob(p))

# Place a panel at an exact position (inches from top-left) with its bold letter.
place_panel <- function(letter, draw_fun, x, top, w, h, fig_h, inset = 0.18) {
  pushViewport(viewport(x = unit(x + inset, "in"), y = unit(fig_h - top, "in"),
                        width = unit(w, "in"), height = unit(h, "in"), just = c("left", "top")))
  draw_fun(); popViewport()
  if (!is.na(letter))
    grid.text(letter, x = unit(x + 0.02, "in"), y = unit(fig_h - top - 0.02, "in"), just = c("left", "top"),
              gp = gpar(fontsize = FS_LETTER, fontface = "bold", fontfamily = FONT))
}
