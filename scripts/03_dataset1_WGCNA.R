# ============================================================================
# 03_dataset1_WGCNA.R
# Weighted gene co-expression network analysis on dataset 1 (Fig 2D, 2E, S2D).
#   input   : VST matrix from 01 (22,863 genes x 18 samples)
#   genes   : top 25% most variable (5,716)
#   network : signed hybrid, TOMType signed, soft power chosen by
#             pickSoftThreshold (RsquaredCut 0.85 -> 18), blockwiseModules with
#             minModuleSize 30, mergeCutHeight 0.25, single block (maxBlockSize 20000)
#   traits  : binary indicators for each genotype and diet; Pearson correlation
#             of module eigengenes with traits, Student p-values (Fig 2D)
#   modules : GO BP and KEGG over-representation per module (organism-wide
#             background, BH p.adjust < 0.05)                        -> Fig 2E
#
# Outputs (results/dataset1_WGCNA/):
#   soft_threshold_selection.csv/.pdf, gene_module_assignments.csv,
#   module_eigengenes.csv, module_trait_correlations.csv, module_trait_pvalues.csv,
#   module_dendrogram.pdf, module_trait_relationships.pdf, module_eigengenes.pdf,
#   Targeted_Module_Enrichment/{GO_BP,KEGG}_<module>_Module.csv (+ plots),
#   Fig2E_KEGG_top5_per_module.csv   (5 pathways per trait-significant module, by Count)
#
# Runtime note: blockwiseModules on 5,716 genes takes a few minutes; the TOM is
# saved to TOM-block.1.RData so re-runs can set loadTOM = TRUE.
# ============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), "..", "R", "common.R"))
need_pkgs("WGCNA", "clusterProfiler", "org.Dr.eg.db", "dplyr", "ggplot2", "enrichplot")
OUT <- RES_D1_WG
allowWGCNAThreads()
set.seed(12345)

# ============================================================================
# 1. LOAD AND FILTER
# ============================================================================
banner("1. Load VST and keep the top 25% most variable genes")
vsd_mat <- read.csv(file.path(RES_D1, "vst_transformed_counts.csv"), row.names = 1, check.names = FALSE)
sample_info <- read.csv(file.path(RES_D1, "sample_metadata.csv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(vsd_mat), sample_info$sample_id))
datExpr <- t(as.matrix(vsd_mat))
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
gene_vars <- apply(datExpr, 2, var)
cutoff <- quantile(gene_vars, probs = 0.75)
datExpr <- datExpr[, gene_vars > cutoff]
kv("genes all / top 25%", c(ncol(vsd_mat), ncol(datExpr)))
check_expected("wgcna_input_genes", ncol(datExpr))

traitData <- data.frame(
  Genotype_EK   = as.numeric(sample_info$genotype == "EK"),
  Genotype_1116 = as.numeric(sample_info$genotype == "pd1116"),
  Genotype_1245 = as.numeric(sample_info$genotype == "pd1245"),
  Diet_HP       = as.numeric(sample_info$diet == "HP"),
  Diet_LP       = as.numeric(sample_info$diet == "LP"))
rownames(traitData) <- sample_info$sample_id

pdf(file.path(OUT, "sample_clustering.pdf"), width = 12, height = 8)
sampleTree <- hclust(dist(datExpr), method = "average")
plotDendroAndColors(sampleTree, numbers2colors(traitData, signed = FALSE), groupLabels = names(traitData),
                    main = "Sample dendrogram and trait heatmap")
invisible(dev.off())

# ============================================================================
# 2. SOFT THRESHOLD
# ============================================================================
banner("2. Soft-thresholding power")
powers <- seq(4, 20, by = 2)
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed hybrid", verbose = 0)
write.csv(sft$fitIndices, file.path(OUT, "soft_threshold_selection.csv"), row.names = FALSE)
print(sft$fitIndices[, c("Power", "SFT.R.sq", "slope", "mean.k.")])
softPower <- sft$powerEstimate
if (is.na(softPower)) {                       # fallback: first power with signed R^2 > 0.8
  r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  softPower <- if (any(r2 > 0.8)) sft$fitIndices$Power[which(r2 > 0.8)[1]] else 16
}
kv("soft power (first R^2 >= 0.85)", softPower)
check_expected("wgcna_soft_power", softPower)
pdf(file.path(OUT, "soft_threshold_selection.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], type = "n", ylim = c(0, 1),
     xlab = "Soft threshold (power)", ylab = "Scale-free topology fit, signed R^2", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, col = "red")
abline(h = c(0.8, 0.85), col = c("red", "blue"), lty = c(1, 2))
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft threshold (power)", ylab = "Mean connectivity", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
invisible(dev.off())

# ============================================================================
# 3. NETWORK AND MODULES
# ============================================================================
banner("3. blockwiseModules (signed hybrid, single block)")
cor <- WGCNA::cor
net <- blockwiseModules(datExpr, power = softPower, networkType = "signed hybrid", TOMType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE,
                        pamRespectsDendro = FALSE, maxBlockSize = 20000,
                        saveTOMs = TRUE, saveTOMFileBase = file.path(OUT, "TOM"), loadTOM = FALSE, verbose = 0)
moduleColors <- labels2colors(net$colors)
modules <- data.frame(gene_id = colnames(datExpr), module_number = net$colors, module_color = moduleColors)
write.csv(modules, file.path(OUT, "gene_module_assignments.csv"), row.names = FALSE)
kv("modules incl. grey", length(unique(moduleColors)))
print(table(moduleColors))
check_expected("wgcna_n_modules_incl_grey", length(unique(moduleColors)))

pdf(file.path(OUT, "module_dendrogram.pdf"), width = 12, height = 8)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
invisible(dev.off())

# ============================================================================
# 4. MODULE EIGENGENES AND TRAIT CORRELATIONS (Fig 2D, S2D)
# ============================================================================
banner("4. Module eigengenes vs traits")
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
write.csv(data.frame(sample_id = rownames(MEs), MEs, check.names = FALSE), file.path(OUT, "module_eigengenes.csv"), row.names = FALSE)
moduleTraitCor    <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
write.csv(moduleTraitCor, file.path(OUT, "module_trait_correlations.csv"))
write.csv(moduleTraitPvalue, file.path(OUT, "module_trait_pvalues.csv"))
non_grey <- rownames(moduleTraitCor) != "MEgrey"
sig_mod  <- rownames(moduleTraitCor)[non_grey][apply(moduleTraitPvalue[non_grey, ] < 0.05, 1, any)]
kv("trait-significant modules (p < 0.05, grey excluded)", sub("^ME", "", sig_mod))
check_expected("wgcna_sig_modules", length(sig_mod))
check_expected("wgcna_star_cells", sum(moduleTraitPvalue[non_grey, ] < 0.05))
print(round(moduleTraitCor, 2))

pdf(file.path(OUT, "module_trait_relationships.pdf"), width = 8, height = 10)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(traitData), yLabels = names(MEs), ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.6, zlim = c(-1, 1), main = "Module-trait relationships")
invisible(dev.off())

pdf(file.path(OUT, "module_eigengenes.pdf"), width = 12, height = 10)
par(mfrow = c(3, 2))
cond <- factor(sample_info$condition, levels = COND_LEVELS)
for (m in names(MEs)) boxplot(MEs[[m]] ~ cond, main = paste("Module", sub("^ME", "", m)), xlab = "", ylab = "Eigengene",
                              col = sub("^ME", "", m), las = 2, outline = FALSE)
invisible(dev.off())

# ============================================================================
# 5. PER-MODULE ENRICHMENT (GO BP, KEGG) -> Fig 2E table
# ============================================================================
banner("5. Per-module GO BP and KEGG (organism-wide background)")
enr_dir <- file.path(OUT, "Targeted_Module_Enrichment"); dir.create(enr_dir, showWarnings = FALSE)
kegg_info <- tryCatch(readLines("https://rest.kegg.jp/info/dre"), error = function(e) "KEGG release query failed")
writeLines(c(paste("Queried:", format(Sys.time())), kegg_info), file.path(OUT, "KEGG_release.txt"))
kegg_tables <- list()
for (color in sort(unique(moduleColors))) {
  genes <- colnames(datExpr)[moduleColors == color]
  write.csv(data.frame(gene_id = genes), file.path(enr_dir, paste0("GeneList_", color, ".csv")), row.names = FALSE)
  map <- suppressMessages(tryCatch(bitr(genes, "SYMBOL", "ENTREZID", org.Dr.eg.db), error = function(e) NULL))
  cat(sprintf("  %-12s %4d genes, %4d with Entrez ID\n", color, length(genes), if (is.null(map)) 0 else nrow(map)))
  if (is.null(map) || nrow(map) == 0) next
  ego <- enrichGO(map$ENTREZID, OrgDb = org.Dr.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(ego) && nrow(as.data.frame(ego))) {
    write.csv(as.data.frame(ego), file.path(enr_dir, paste0("GO_BP_", color, "_Module.csv")), row.names = FALSE)
    pdf(file.path(enr_dir, paste0("GO_Plot_", color, "_Module.pdf")), width = 10, height = 8)
    print(dotplot(ego, showCategory = 20, title = paste("GO BP:", color, "module"))); invisible(dev.off())
  }
  kk <- enrichKEGG(map$ENTREZID, organism = "dre", pvalueCutoff = 0.05)
  if (!is.null(kk) && nrow(as.data.frame(kk))) {
    kk <- setReadable(kk, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")
    d <- as.data.frame(kk); kegg_tables[[color]] <- d
    write.csv(d, file.path(enr_dir, paste0("KEGG_", color, "_Module.csv")), row.names = FALSE)
    pdf(file.path(enr_dir, paste0("KEGG_Plot_", color, "_Module.pdf")), width = 10, height = 8)
    print(dotplot(kk, showCategory = 20, title = paste("KEGG:", color, "module"))); invisible(dev.off())
  }
}

# Fig 2E: five pathways per trait-significant module, chosen by gene Count among
# significant pathways (as published); grey excluded.
fig2e_modules <- sub("^ME", "", sig_mod)
fig2e <- bind_rows(lapply(fig2e_modules, function(m) {
  d <- kegg_tables[[m]]; if (is.null(d) || !nrow(d)) return(NULL)
  d <- d[order(-d$Count), ][seq_len(min(5, nrow(d))), ]; d$Cluster <- m; d }))
fig2e$GeneRatioNum <- sapply(strsplit(fig2e$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
write.csv(fig2e, file.path(OUT, "Fig2E_KEGG_top5_per_module.csv"), row.names = FALSE)
kv("Fig 2E rows (modules x 5)", nrow(fig2e))
check_expected("fig2E_dots", nrow(fig2e))

banner("6. Output verification")
verify_files(file.path(OUT, c("soft_threshold_selection.csv", "gene_module_assignments.csv", "module_eigengenes.csv",
                              "module_trait_correlations.csv", "module_trait_pvalues.csv", "Fig2E_KEGG_top5_per_module.csv",
                              "module_dendrogram.pdf", "module_trait_relationships.pdf", "module_eigengenes.pdf")))
write_session()
cat("\nDone: 03_dataset1_WGCNA.R\n")
