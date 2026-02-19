# ============================================================================
# FINAL PUBLICATION SCRIPT: Enrichment Analysis with Fixed X-Axis
# ============================================================================

# 1. SETUP
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(org.Dr.eg.db)
library(ggplot2)
library(stringr)
library(pheatmap)
library(tibble)
library(RColorBrewer)

output_dir <- "../analysis_results"

# ============================================================================
# 2. ROBUST FUNCTIONS
# ============================================================================

# A. Load Genes
get_genes <- function(filename, direction) {
  fpath <- file.path(output_dir, filename)
  if(!file.exists(fpath)) { message("Warning: File not found -> ", filename); return(character(0)) }
  df <- read.csv(fpath)
  if(direction == "up") return(df$gene_id[df$log2FoldChange > 0])
  else return(df$gene_id[df$log2FoldChange < 0])
}

# B. Robust ID Conversion (Symbol -> Entrez)
to_entrez <- function(genes) {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  if(length(genes) == 0) return(NULL)
  
  # Try Symbol
  map1 <- suppressMessages(tryCatch(bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dr.eg.db), error=function(e) NULL))
  # Try Ensembl
  map2 <- suppressMessages(tryCatch(bitr(genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Dr.eg.db), error=function(e) NULL))
  
  ids <- unique(c(map1$ENTREZID, map2$ENTREZID))
  return(ids)
}

# C. Clean Symbols (Remove Ensembl IDs)
clean_symbols <- function(gene_list) {
  return(gene_list[!grepl("^ENSDARG", gene_list)])
}

# D. CUSTOM PLOTTING FUNCTION (With Text Wrapping & Larger Font)
custom_dotplot <- function(enrich_result, title, expected_levels) {
  # 1. Extract Data
  if (is.null(enrich_result)) { return(NULL) }
  df <- as.data.frame(enrich_result)
  
  if (nrow(df) == 0) {
    message("No significant terms found for: ", title)
    return(NULL)
  }
  
  # 2. Select Top 10 terms per group
  top_terms <- df %>% 
    group_by(Cluster) %>% 
    slice_min(order_by = p.adjust, n = 10) %>% 
    pull(Description)
  
  plot_data <- df %>% filter(Description %in% top_terms)
  
  # 3. FORCE FACTOR LEVELS
  plot_data$Cluster <- factor(plot_data$Cluster, levels = expected_levels)
  
  # 4. GGPLOT with Text Wrapping
  p <- ggplot(plot_data, aes(x = Cluster, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    
    # --- FIX 1: FORCE ALL COLUMNS TO SHOW ---
    scale_x_discrete(drop = FALSE) +  
    
    # --- FIX 2: WRAP LONG LABELS ---
    # width = 50 means "break line after ~50 characters"
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) + 
    
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL, color = "FDR", size = "Count") +
    
    # --- FIX 3: INCREASE FONT SIZE ---
    theme(
      # X-axis: Angle 45, Bold, Size 11
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 11),
      
      # Y-axis: Size 12 (Larger), Black color
      axis.text.y = element_text(size = 12, color = "black"),
      
      # Title: Large and Bold
      plot.title = element_text(size = 16, face = "bold"),
      
      # Add a little padding to the right of Y labels so they don't touch the dots
      axis.title.y = element_blank(), 
      panel.grid.major = element_line(color = "grey90")
    )
  
  return(p)
}

# ============================================================================
# 3. PREPARE DATA
# ============================================================================
cat("Loading Background Universe...\n")
counts <- read.csv("../analysis_results/normalized_counts.csv", row.names = 1)
univ_sym <- rownames(counts)
univ_ent <- to_entrez(univ_sym)

# Load Lists (HP)
cat("Loading HP Lists...\n")
hp_sym_1116_up   <- get_genes("DESeq2_significant_HP_1116_vs_EK.csv", "up")
hp_sym_1116_down <- get_genes("DESeq2_significant_HP_1116_vs_EK.csv", "down")
hp_sym_1245_up   <- get_genes("DESeq2_significant_HP_1245_vs_EK.csv", "up")
hp_sym_1245_down <- get_genes("DESeq2_significant_HP_1245_vs_EK.csv", "down")

# Load Lists (LP)
cat("Loading LP Lists...\n")
lp_sym_1116_up   <- get_genes("DESeq2_significant_LP_1116_vs_EK.csv", "up")
lp_sym_1116_down <- get_genes("DESeq2_significant_LP_1116_vs_EK.csv", "down")
lp_sym_1245_up   <- get_genes("DESeq2_significant_LP_1245_vs_EK.csv", "up")
lp_sym_1245_down <- get_genes("DESeq2_significant_LP_1245_vs_EK.csv", "down")

# Load Lists (HP vs LP)
cat("Loading HP vs LP Lists...\n")
hpvslp_sym_EK_up   <- get_genes("DESeq2_significant_EK_HP_vs_LP.csv", "up")
hpvslp_sym_EK_down <- get_genes("DESeq2_significant_EK_HP_vs_LP.csv", "down")
hpvslp_sym_1116_up   <- get_genes("DESeq2_significant_1116_HP_vs_LP.csv", "up")
hpvslp_sym_1116_down <- get_genes("DESeq2_significant_1116_HP_vs_LP.csv", "down")
hpvslp_sym_1245_up   <- get_genes("DESeq2_significant_1245_HP_vs_LP.csv", "up")
hpvslp_sym_1245_down <- get_genes("DESeq2_significant_1245_HP_vs_LP.csv", "down")

# Load Lists (HP vs LP)
cat("Loading HP vs LP Lists...\n")
hp_1116_unique_up   <- get_genes("genes_unique_to_1116_HP_table.csv", "up")
hp_1116_unique_down <- get_genes("genes_unique_to_1116_HP_table.csv", "down")
hp_1245_unique_up   <- get_genes("genes_unique_to_1245_HP_table.csv", "up")
hp_1245_unique_down <- get_genes("genes_unique_to_1245_HP_table.csv", "down")
lp_1116_unique_up   <- get_genes("genes_unique_to_1116_LP_table.csv", "up")
lp_1116_unique_down <- get_genes("genes_unique_to_1116_LP_table.csv", "down")
lp_1245_unique_up   <- get_genes("genes_unique_to_1245_LP_table.csv", "up")
lp_1245_unique_down <- get_genes("genes_unique_to_1245_LP_table.csv", "down")

# Define Factor Levels (Order of columns)
levels_HP <- c("1116_vs_EK_HP_Up", "1116_vs_EK_HP_Down", "1245_vs_EK_HP_Up", "1245_vs_EK_HP_Down")
levels_LP <- c("1116_vs_EK_LP_Up", "1116_vs_EK_LP_Down", "1245_vs_EK_LP_Up", "1245_vs_EK_LP_Down")
levels_HPvsLP <- c("EK_HPvsLP_Up", "EK_HPvsLP_Down", 
                   "1116_HPvsLP_Up", "1116_HPvsLP_Down", 
                   "1245_HPvsLP_Up", "1245_HPvsLP_Down")
levels_unique <- c("HP_1116_unique_Up", "HP_1116_unique_Down",
                   "HP_1245_unique_Up", "HP_1245_unique_Down",
                   "LP_1116_unique_Up", "LP_1116_unique_Down",
                   "LP_1245_unique_Up", "LP_1245_unique_Down")

levels_unique_1245 <- c("HP_1245_unique_Up", "HP_1245_unique_Down", "LP_1245_unique_Up", "LP_1245_unique_Down")


# ============================================================================
# 4. RUN ANALYSIS (KEGG)
# ============================================================================
cat("\nRunning KEGG...\n")

# HP
kegg_input_HP <- list(
  "1116_vs_EK_HP_Up"   = to_entrez(hp_sym_1116_up),
  "1116_vs_EK_HP_Down" = to_entrez(hp_sym_1116_down),
  "1245_vs_EK_HP_Up"   = to_entrez(hp_sym_1245_up),
  "1245_vs_EK_HP_Down" = to_entrez(hp_sym_1245_down)
)
ck_kegg_HP <- compareCluster(kegg_input_HP, fun="enrichKEGG", organism="dre", universe=univ_ent, pvalueCutoff=0.05)
ck_kegg_HP <- setReadable(ck_kegg_HP, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")

# LP
kegg_input_LP <- list(
  "1116_vs_EK_LP_Up"   = to_entrez(lp_sym_1116_up),
  "1116_vs_EK_LP_Down" = to_entrez(lp_sym_1116_down),
  "1245_vs_EK_LP_Up"   = to_entrez(lp_sym_1245_up),
  "1245_vs_EK_LP_Down" = to_entrez(lp_sym_1245_down)
)
ck_kegg_LP <- compareCluster(kegg_input_LP, fun="enrichKEGG", organism="dre", universe=univ_ent, pvalueCutoff=0.05)
ck_kegg_LP <- setReadable(ck_kegg_LP, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")

# HP vs LP
kegg_input_HPvsLP <- list(
  "EK_HPvsLP_Up"   = to_entrez(hpvslp_sym_EK_up),
  "EK_HPvsLP_Down" = to_entrez(hpvslp_sym_EK_down),
  "1116_HPvsLP_Up"   = to_entrez(hpvslp_sym_1116_up),
  "1116_HPvsLP_Down" = to_entrez(hpvslp_sym_1116_down),
  "1245_HPvsLP_Up"   = to_entrez(hpvslp_sym_1245_up),
  "1245_HPvsLP_Down" = to_entrez(hpvslp_sym_1245_down)
)
ck_kegg_HPvsLP <- compareCluster(kegg_input_HPvsLP, fun="enrichKEGG", organism="dre", universe=univ_ent, pvalueCutoff=0.05)
ck_kegg_HPvsLP <- setReadable(ck_kegg_HPvsLP, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")

# 1245 unique HP and LP --> only KEGG ran for now
kegg_input_1245unique <- list(
  "HP_1245_unique_Up"   = to_entrez(hp_1245_unique_up),
  "HP_1245_unique_Down" = to_entrez(hp_1245_unique_down),
  "LP_1245_unique_Up"   = to_entrez(lp_1245_unique_up),
  "LP_1245_unique_Down" = to_entrez(lp_1245_unique_down)
)
ck_kegg_1245unique <- compareCluster(kegg_input_1245unique, fun="enrichKEGG", organism="dre", universe=univ_ent, pvalueCutoff=0.05)
ck_kegg_1245unique <- setReadable(ck_kegg_1245unique, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")


# ============================================================================
# 5. RUN ANALYSIS (GO BP)
# ============================================================================
cat("Running GO BP...\n")

# HP
go_input_HP <- list(
  "1116_vs_EK_HP_Up"   = clean_symbols(hp_sym_1116_up),
  "1116_vs_EK_HP_Down" = clean_symbols(hp_sym_1116_down),
  "1245_vs_EK_HP_Up"   = clean_symbols(hp_sym_1245_up),
  "1245_vs_EK_HP_Down" = clean_symbols(hp_sym_1245_down)
)
ck_go_HP <- compareCluster(go_input_HP, fun="enrichGO", OrgDb=org.Dr.eg.db, keyType="SYMBOL", ont="BP", universe=univ_sym, pvalueCutoff=0.05)
ck_go_HP_simp <- simplify(ck_go_HP, cutoff=0.6, by="p.adjust", select_fun=min)

# LP
go_input_LP <- list(
  "1116_vs_EK_LP_Up"   = clean_symbols(lp_sym_1116_up),
  "1116_vs_EK_LP_Down" = clean_symbols(lp_sym_1116_down),
  "1245_vs_EK_LP_Up"   = clean_symbols(lp_sym_1245_up),
  "1245_vs_EK_LP_Down" = clean_symbols(lp_sym_1245_down)
)
ck_go_LP <- compareCluster(go_input_LP, fun="enrichGO", OrgDb=org.Dr.eg.db, keyType="SYMBOL", ont="BP", universe=univ_sym, pvalueCutoff=0.05)
ck_go_LP_simp <- simplify(ck_go_LP, cutoff=0.6, by="p.adjust", select_fun=min)

# HPvsLP
go_input_HPvsLP <- list(
  "EK_HPvsLP_Up"   = clean_symbols(hpvslp_sym_EK_up),
  "EK_HPvsLP_Down" = clean_symbols(hpvslp_sym_EK_down),
  "1116_HPvsLP_Up"   = clean_symbols(hpvslp_sym_1116_up),
  "1116_HPvsLP_Down" = clean_symbols(hpvslp_sym_1116_down),
  "1245_HPvsLP_Up"   = clean_symbols(hpvslp_sym_1245_up),
  "1245_HPvsLP_Down" = clean_symbols(hpvslp_sym_1245_down)
)
ck_go_HPvsLP <- compareCluster(go_input_HPvsLP, fun="enrichGO", OrgDb=org.Dr.eg.db, keyType="SYMBOL", ont="BP", universe=univ_sym, pvalueCutoff=0.05)
ck_go_HPvsLP_simp <- simplify(ck_go_HPvsLP, cutoff=0.6, by="p.adjust", select_fun=min)

# ============================================================================
# 6. GENERATE PLOTS
# ============================================================================
cat("\nSaving Plots...\n")

pdf(file.path(output_dir, "Final_Enrichment_Plots.pdf"), width = 12, height = 18)

# Use our CUSTOM function
print(custom_dotplot(ck_kegg_HP, "KEGG Pathways: High Protein", levels_HP))
print(custom_dotplot(ck_kegg_LP, "KEGG Pathways: Low Protein", levels_LP))
print(custom_dotplot(ck_kegg_HPvsLP, "KEGG Pathways: High vs Low Protein", levels_HPvsLP))
print(custom_dotplot(ck_go_HP_simp, "GO BP: High Protein", levels_HP))
print(custom_dotplot(ck_go_LP_simp, "GO BP: Low Protein", levels_LP))
print(custom_dotplot(ck_go_HPvsLP_simp, "GO BP: High vs Low Protein", levels_HPvsLP))

dev.off()
cat("Done! Check 'Final_Enrichment_Plots.pdf'\n")

# ============================================================================
# 7. SAVE ENRICHMENT RESULTS TO CSV/XLSX
# ============================================================================
cat("Exporting Enrichment Results to CSV...\n")

# Helper function to save if object exists
save_enrich <- function(enrich_obj, filename) {
  if (!is.null(enrich_obj)) {
    # Convert to standard dataframe
    df <- as.data.frame(enrich_obj)
    write.csv(df, file.path(output_dir, filename), row.names = FALSE)
    cat("Saved:", filename, "\n")
  }
}

# Save the 4 main comparisons
save_enrich(ck_kegg_HP, "Enrichment_KEGG_HP.csv")
save_enrich(ck_kegg_LP, "Enrichment_KEGG_LP.csv")
save_enrich(ck_kegg_HPvsLP, "Enrichment_KEGG_HPvsLP.csv")
save_enrich(ck_go_HP_simp, "Enrichment_GO_HP.csv")
save_enrich(ck_go_LP_simp, "Enrichment_GO_LP.csv")
save_enrich(ck_go_HPvsLP_simp, "Enrichment_GO_HPvsLP.csv")
