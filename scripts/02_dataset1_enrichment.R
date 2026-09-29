# ============================================================================
# 02_dataset1_enrichment.R
# Over-representation analysis of dataset-1 DEG lists (clusterProfiler).
#   GO BP  : enrichGO, keyType SYMBOL, background = all 22,863 expressed genes,
#            BH p.adjust < 0.05, then simplify(cutoff 0.6)          -> Fig 5A / 5B
#   KEGG   : enrichKEGG organism "dre", background = expressed genes mapped to
#            Entrez, BH p.adjust < 0.05                              -> Fig S3A
# Gene sets: up / down DEGs for pd1116 vs EK and pd1245 vs EK (per diet), HP vs
# LP within genotype, and DEGs unique to pd1245 vs EK (not in pd1116 vs EK).
#
# NOTE ON KEGG: enrichKEGG downloads the current KEGG release at run time, so
# pathway membership can drift between runs. The release used is recorded in
# results/dataset1_enrichment/KEGG_release.txt. The panel published as Fig S3A
# was produced from the Feb 2026 release (11 pathways).
#
# Outputs (results/dataset1_enrichment/):
#   Enrichment_GO_{HP,LP,HPvsLP}.csv        (simplified compareCluster tables)
#   Enrichment_KEGG_{HP,LP,HPvsLP,1245unique}.csv
#   compareCluster_*.rds                     (objects, for the figure scripts)
#   Enrichment_dotplots.pdf                  (quick-look; publication: 06_figure5AB.R)
# ============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), "..", "R", "common.R"))
need_pkgs("clusterProfiler", "org.Dr.eg.db", "dplyr", "ggplot2", "stringr")
OUT <- RES_D1_ENR

# ============================================================================
# 1. GENE LISTS
# ============================================================================
banner("1. Load DEG lists and background")
get_genes <- function(file, direction) {
  d <- read.csv(file.path(RES_D1, file), stringsAsFactors = FALSE)
  if (direction == "up") d$gene_id[d$log2FoldChange > 0] else d$gene_id[d$log2FoldChange < 0]
}
univ_sym <- rownames(read.csv(file.path(RES_D1, "normalized_counts.csv"), row.names = 1))
kv("background (expressed genes, symbols)", length(univ_sym))

# Symbol -> Entrez, with an Ensembl fallback for genes kept under their Ensembl ID
to_entrez <- function(genes) {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  if (!length(genes)) return(NULL)
  m1 <- suppressMessages(tryCatch(bitr(genes, "SYMBOL", "ENTREZID", org.Dr.eg.db), error = function(e) NULL))
  m2 <- suppressMessages(tryCatch(bitr(genes, "ENSEMBL", "ENTREZID", org.Dr.eg.db), error = function(e) NULL))
  unique(c(m1$ENTREZID, m2$ENTREZID))
}
clean_symbols <- function(x) x[!grepl("^ENSDARG", x)]
univ_ent <- to_entrez(univ_sym)
kv("background mapped to Entrez (KEGG universe)", length(univ_ent))
kv("org.Dr.eg.db version", as.character(packageVersion("org.Dr.eg.db")))

lists <- list(
  HP = list("1116_vs_EK_HP_Up"   = get_genes("DESeq2_significant_HP_1116_vs_EK.csv", "up"),
            "1116_vs_EK_HP_Down" = get_genes("DESeq2_significant_HP_1116_vs_EK.csv", "down"),
            "1245_vs_EK_HP_Up"   = get_genes("DESeq2_significant_HP_1245_vs_EK.csv", "up"),
            "1245_vs_EK_HP_Down" = get_genes("DESeq2_significant_HP_1245_vs_EK.csv", "down")),
  LP = list("1116_vs_EK_LP_Up"   = get_genes("DESeq2_significant_LP_1116_vs_EK.csv", "up"),
            "1116_vs_EK_LP_Down" = get_genes("DESeq2_significant_LP_1116_vs_EK.csv", "down"),
            "1245_vs_EK_LP_Up"   = get_genes("DESeq2_significant_LP_1245_vs_EK.csv", "up"),
            "1245_vs_EK_LP_Down" = get_genes("DESeq2_significant_LP_1245_vs_EK.csv", "down")),
  HPvsLP = list("EK_HPvsLP_Up"     = get_genes("DESeq2_significant_EK_HP_vs_LP.csv", "up"),
                "EK_HPvsLP_Down"   = get_genes("DESeq2_significant_EK_HP_vs_LP.csv", "down"),
                "1116_HPvsLP_Up"   = get_genes("DESeq2_significant_1116_HP_vs_LP.csv", "up"),
                "1116_HPvsLP_Down" = get_genes("DESeq2_significant_1116_HP_vs_LP.csv", "down"),
                "1245_HPvsLP_Up"   = get_genes("DESeq2_significant_1245_HP_vs_LP.csv", "up"),
                "1245_HPvsLP_Down" = get_genes("DESeq2_significant_1245_HP_vs_LP.csv", "down")),
  `1245unique` = list("HP_1245_unique_Up"   = get_genes("genes_unique_to_1245_HP_table.csv", "up"),
                      "HP_1245_unique_Down" = get_genes("genes_unique_to_1245_HP_table.csv", "down"),
                      "LP_1245_unique_Up"   = get_genes("genes_unique_to_1245_LP_table.csv", "up"),
                      "LP_1245_unique_Down" = get_genes("genes_unique_to_1245_LP_table.csv", "down")))
for (set in names(lists)) { cat("  ", set, ":\n"); for (k in names(lists[[set]])) kv(paste0("    ", k), length(lists[[set]][[k]])) }

# ============================================================================
# 2. KEGG (all four sets)
# ============================================================================
banner("2. KEGG over-representation (enrichKEGG, dre)")
kegg_info <- tryCatch(readLines("https://rest.kegg.jp/info/dre"), error = function(e) "KEGG release query failed")
writeLines(c(paste("Queried:", format(Sys.time())), kegg_info), file.path(OUT, "KEGG_release.txt"))
kv("KEGG release", kegg_info[1])
ck_kegg <- list()
for (set in names(lists)) {
  ent <- lapply(lists[[set]], to_entrez)
  ck  <- compareCluster(ent, fun = "enrichKEGG", organism = "dre", universe = univ_ent, pvalueCutoff = 0.05)
  ck  <- setReadable(ck, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")
  ck_kegg[[set]] <- ck
  df <- as.data.frame(ck)
  write.csv(df, file.path(OUT, paste0("Enrichment_KEGG_", set, ".csv")), row.names = FALSE)
  saveRDS(ck, file.path(OUT, paste0("compareCluster_KEGG_", set, ".rds")))
  kv(paste("KEGG", set, "significant rows"), nrow(df))
  if (nrow(df)) print(table(df$Cluster))
}

# ============================================================================
# 3. GO BP (HP, LP, HPvsLP) + simplify
# ============================================================================
banner("3. GO BP over-representation (enrichGO, SYMBOL) + simplify(0.6)")
ck_go <- list()
for (set in c("HP", "LP", "HPvsLP")) {
  sy <- lapply(lists[[set]], clean_symbols)
  ck <- compareCluster(sy, fun = "enrichGO", OrgDb = org.Dr.eg.db, keyType = "SYMBOL", ont = "BP",
                       universe = univ_sym, pvalueCutoff = 0.05)
  cks <- simplify(ck, cutoff = 0.6, by = "p.adjust", select_fun = min)
  ck_go[[set]] <- cks
  df <- as.data.frame(cks)
  write.csv(df, file.path(OUT, paste0("Enrichment_GO_", set, ".csv")), row.names = FALSE)
  saveRDS(cks, file.path(OUT, paste0("compareCluster_GO_", set, ".rds")))
  kv(paste("GO BP", set, "rows before / after simplify"), c(nrow(as.data.frame(ck)), nrow(df)))
  if (nrow(df)) print(table(df$Cluster))
}

# Terms that the Fig 5A/5B rule (top 10 per cluster by p.adjust) would draw
top_terms <- function(df) { if (!nrow(df)) return(character(0))
  df %>% group_by(Cluster) %>% slice_min(order_by = p.adjust, n = 10) %>% pull(Description) %>% unique() }
check_expected("fig5A_terms", length(top_terms(as.data.frame(ck_go$HP))))
check_expected("fig5B_terms", length(top_terms(as.data.frame(ck_go$LP))))

# ============================================================================
# 4. QUICK-LOOK PLOTS
# ============================================================================
banner("4. Quick-look dot plots")
quick_dot <- function(ck, title, levels) {
  df <- as.data.frame(ck); if (!nrow(df)) return(NULL)
  keep <- df %>% group_by(Cluster) %>% slice_min(order_by = p.adjust, n = 10) %>% pull(Description)
  pd <- df %>% filter(Description %in% keep); pd$Cluster <- factor(pd$Cluster, levels = levels)
  ggplot(pd, aes(Cluster, Description, size = Count, color = p.adjust)) + geom_point() +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(labels = function(x) str_wrap(x, 50)) +
    scale_color_gradient(low = "red", high = "blue") + theme_bw() +
    labs(title = title, x = NULL, y = NULL, color = "FDR") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
pdf(file.path(OUT, "Enrichment_dotplots.pdf"), width = 11, height = 14)
for (set in names(lists)) print(quick_dot(ck_kegg[[set]], paste("KEGG:", set), names(lists[[set]])))
for (set in names(ck_go)) print(quick_dot(ck_go[[set]], paste("GO BP:", set), names(lists[[set]])))
invisible(dev.off())

banner("5. Output verification")
verify_files(file.path(OUT, c(paste0("Enrichment_KEGG_", names(lists), ".csv"),
                              paste0("Enrichment_GO_", names(ck_go), ".csv"),
                              "Enrichment_dotplots.pdf", "KEGG_release.txt")))
write_session()
cat("\nDone: 02_dataset1_enrichment.R\n")
