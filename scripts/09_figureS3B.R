# ============================================================================
# 09_figureS3B.R  --  Figure S3B: lipid / cholesterol / vitamin-absorption genes
# 25 curated genes (curated_genes/figS3B_lipid_vitamin.tsv) x 18 samples,
# per-gene z-score of VST across all samples; samples in metadata order
# (genotype-major, diet-minor), annotated by genotype and diet.
# (Fig S3A, the KEGG dot plot of pd1245-unique DEGs, is drawn by 02's quick-look
#  PDF from Enrichment_KEGG_1245unique.csv; see README on KEGG-release drift.)
# ============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), "..", "R", "common.R"))
source(file.path(ROOT, "R", "fig_style.R"))
need_pkgs("ComplexHeatmap", "circlize", "RColorBrewer")
ht_opt(message = FALSE)
OUT <- file.path(FIGURES, "FigureS3"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIG_W <- 7.5; FIG_H <- 4.2

cur <- read.delim(file.path(CURATED, "figS3B_lipid_vitamin.tsv"), stringsAsFactors = FALSE)
cur$group <- factor(cur$group, levels = unique(cur$group))
vst <- as.matrix(read.csv(file.path(RES_D1, "vst_transformed_counts.csv"), row.names = 1, check.names = FALSE))
annot <- read.csv(file.path(RES_D1, "gene_annotation.csv"), stringsAsFactors = FALSE)
sample_info <- read.csv(file.path(RES_D1, "sample_metadata.csv"), stringsAsFactors = FALSE)
row <- annot$symbol_used[match(cur$ensembl_gene_id, annot$ensembl_id)]
row[is.na(row)] <- annot$symbol_used[match(tolower(cur$query_symbol[is.na(row)]), tolower(annot$symbol_used))]
stopifnot(!any(is.na(row)), all(row %in% rownames(vst)))
z <- t(scale(t(vst[row, sample_info$sample_id]))); rownames(z) <- cur$query_symbol
banner("S3B"); kv("genes x samples", paste(dim(z), collapse = " x ")); print(table(cur$group)); kv("z range", sprintf("%.2f", range(z)))
check_expected("figS3B_genes", nrow(z))
write.csv(data.frame(group = cur$group, gene = cur$query_symbol, ensembl_gene_id = cur$ensembl_gene_id, z, check.names = FALSE),
          file.path(RES_HM, "figS3B_zscore_persample_matrix.csv"), row.names = FALSE)

cat_colors <- setNames(RColorBrewer::brewer.pal(3, "Set2"), levels(cur$group))
lg <- function(t) list(title_gp = gpar(fontsize = FS_SMALL, fontface = "bold", fontfamily = FONT),
                       labels_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT), grid_width = unit(2.5, "mm"), grid_height = unit(2.5, "mm"))
ht <- Heatmap(z, name = "Z-score", col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT),
              show_column_names = FALSE, column_title = "Lipid/Cholesterol Metabolism + Vitamin Digestion",
              column_title_gp = gpar(fontsize = FS_MAIN, fontface = "bold", fontfamily = FONT),
              row_split = cur$group, row_title = NULL, row_gap = unit(0.8, "mm"), rect_gp = gpar(col = NA),
              top_annotation = HeatmapAnnotation(Genotype = sample_info$genotype, Diet = sample_info$diet,
                                                 col = list(Genotype = GENO_COLORS, Diet = DIET_COLORS), simple_anno_size = unit(2.5, "mm"),
                                                 annotation_name_side = "left", annotation_name_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT),
                                                 annotation_legend_param = list(Genotype = lg(), Diet = lg())),
              right_annotation = rowAnnotation(Category = cur$group, col = list(Category = cat_colors), simple_anno_size = unit(2.5, "mm"),
                                               show_annotation_name = FALSE, annotation_legend_param = list(Category = lg())),
              heatmap_legend_param = list(title = "Expression\n(Z-score)", at = c(-2, -1, 0, 1, 2),
                                          title_gp = gpar(fontsize = FS_SMALL, fontface = "bold", fontfamily = FONT),
                                          labels_gp = gpar(fontsize = FS_SMALL, fontfamily = FONT), grid_width = unit(2.5, "mm"), legend_height = unit(1.1, "in")))
save_both(function() { pushViewport(viewport(gp = gpar(fontfamily = FONT)))
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE, newpage = FALSE, padding = unit(c(1, 1, 1, 1), "mm"))
  popViewport() }, "FigS3B_lipid_heatmap", FIG_W, FIG_H, OUT)
verify_files(file.path(OUT, c("FigS3B_lipid_heatmap.pdf", "FigS3B_lipid_heatmap.tiff")))
write_session()
cat("\nDone: 09_figureS3B.R\n")
