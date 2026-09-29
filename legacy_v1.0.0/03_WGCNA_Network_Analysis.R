# ============================================================================
# WGCNA Analysis for pllp RNA-seq Data
# Weighted Gene Co-expression Network Analysis
# ============================================================================

library(WGCNA)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Dr.eg.db)
library(ggplot2)
library(dplyr)
library(enrichplot)
library(tidyverse)
library(stringr)


# Allow multi-threading
allowWGCNAThreads()

# Set working directory
output_dir <- "../analysis_results/WGCNA"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Load normalized counts from previous analysis
vsd_mat <- read.csv("../analysis_results/vst_transformed_counts.csv", row.names = 1)
sample_info <- read.csv("../analysis_results/sample_metadata.csv")

# Prepare expression data for WGCNA (samples as rows, genes as columns)
datExpr <- t(vsd_mat)
rownames(datExpr) <- colnames(vsd_mat)

# Check for genes and samples with too many missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# ============================================================================
# 1b. FILTER FOR TOP 25% MOST VARIABLE GENES (HVGs)
# ============================================================================

cat("Filtering for top 25% most variable genes...\n")

# 1. Calculate variance for each gene (columns)
gene_vars <- apply(datExpr, 2, var)

# 2. Determine the percentile cutoff (Top 25% = 75th percentile)
#    (Use '0.5' for top 50% if you want to be less strict)
cutoff <- quantile(gene_vars, probs = 0.75) 

# 3. Subset the data
datExpr_all <- datExpr # Keep backup of full data if needed
datExpr <- datExpr[, gene_vars > cutoff]

cat("Original genes:", ncol(datExpr_all), "\n") #Original genes: 22863
cat("Filtered genes (Top 25%):", ncol(datExpr), "\n") #Filtered genes (Top 25%): 5716 

# ============================================================================
# 2. SAMPLE CLUSTERING TO DETECT OUTLIERS
# ============================================================================

pdf(file.path(output_dir, "sample_clustering.pdf"), width = 12, height = 9)

sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot trait heatmap
traitData <- data.frame(
  Genotype_EK = as.numeric(sample_info$genotype == "EK"),
  Genotype_1116 = as.numeric(sample_info$genotype == "pd1116"),
  Genotype_1245 = as.numeric(sample_info$genotype == "pd1245"),
  Diet_HP = as.numeric(sample_info$diet == "HP"),
  Diet_LP = as.numeric(sample_info$diet == "LP")
)

traitColors <- numbers2colors(traitData, signed = FALSE)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(traitData),
                    main = "Sample dendrogram and trait heatmap")

dev.off()

# ============================================================================
# 3. CHOOSE SOFT-THRESHOLDING POWER
# ============================================================================

# 1. Define powers explicitly first (so you can reuse the variable)
powers <- seq(4, 20, by = 2)

# 2. Run the analysis
sft <- pickSoftThreshold(datExpr, 
                         powerVector = powers, 
                         networkType = "signed hybrid", 
                         verbose = 5)

# 3. Plot Results
pdf(file.path(output_dir, "soft_threshold_selection.pdf"), width = 12, height = 5)

par(mfrow = c(1, 2))
cex1 <- 0.9

# PLOT 1: Scale Independence
# Use 'ylim' to make sure you see the top of the curve
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", 
     main = "Scale independence",
     ylim = c(0, 1)) 

text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")

# Add reference lines (Standard 0.8, Strict 0.9)
abline(h = 0.80, col = "red")
abline(h = 0.90, col = "blue", lty = 2)

# PLOT 2: Mean Connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n",
     main = "Mean connectivity")

text(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     labels = powers, cex = cex1, col = "red")

dev.off()

# 4. Smart Power Selection Logic
# Sometimes the auto-picker returns NA if it doesn't hit 0.9 exactly.
# This logic picks the first power > 0.8 manually if the auto-picker fails.
softPower <- sft$powerEstimate
#WGCNA Selected soft power: 18 
if (is.na(softPower)) {
  # Calculate R^2 manually
  r2_values <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  
  # Find first index > 0.8
  pass_threshold <- which(r2_values > 0.8)
  
  if (length(pass_threshold) > 0) {
    softPower <- sft$fitIndices[pass_threshold[1], "Power"]
    cat("Manual Selection: Power", softPower, "is the first to cross R^2 > 0.8.\n")
  } else {
    softPower <- 16 # Fallback based on your previous data
    cat("Warning: No power reached 0.8. Defaulting to 16.\n")
  }
} else {
  cat("WGCNA Selected soft power:", softPower, "\n")
}

# Save power selection info
write.csv(sft$fitIndices, file.path(output_dir, "soft_threshold_selection.csv"), row.names = FALSE)

# ============================================================================
# 4. CONSTRUCT CO-EXPRESSION NETWORK
# ============================================================================

# One-step network construction
cat("Constructing network...\n")
# 1. FIX: Force R to use WGCNA's correlation function
cor <- WGCNA::cor

# 2. RUN: Add 'loadTOM = TRUE' to skip re-calculating the blocks
net <- blockwiseModules(datExpr, 
                        # 1. Critical Biological Parameters (KEEP THESE)
                        power = softPower,
                        networkType = "signed hybrid", 
                        TOMType = "signed", 
                        minModuleSize = 30,
                        mergeCutHeight = 0.25,      # Controls merging of similar modules
                        
                        # 2. Output Format (KEEP THIS for your script to work)
                        numericLabels = TRUE,       # Your script uses labels2colors() later, so keep this TRUE
                        pamRespectsDendro = FALSE,  # Standard setting
                        # --- THE CHANGE: Set this to TRUE to reload ---
                        loadTOM = FALSE,             # <--- CHANGE THIS TO TRUE when re-running to save time
                        # ----------------------------------------------
                        # 3. The "Single Block" Fix
                        maxBlockSize = 20000,       # Set this > 5716 (your gene count)
                        
                        # 4. Parameters you can OMIT or Simplfy
                        # reassignThreshold: Not needed for single block
                        # loadTOM: REMOVED (caused your error)
                        
                        # Optional: Keep saving if you want a backup, but not strictly necessary for small data
                        saveTOMs = TRUE,
                        saveTOMFileBase = file.path(output_dir, "TOM"),
                        verbose = 3)

# Save module assignments
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
modules <- data.frame(
  gene_id = colnames(datExpr),
  module_number = moduleLabels,
  module_color = moduleColors
)
write.csv(modules, file.path(output_dir, "gene_module_assignments.csv"), row.names = FALSE)

cat("Number of modules:", length(unique(moduleColors)), "\n")
cat("Module sizes:\n")
print(table(moduleColors))

# ============================================================================
# 5. VISUALIZE THE NETWORK
# ============================================================================

pdf(file.path(output_dir, "module_dendrogram.pdf"), width = 12, height = 9)

# Plot the dendrogram and the module colors
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

dev.off()

# ============================================================================
# 6. RELATE MODULES TO EXTERNAL TRAITS
# ============================================================================

# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Correlate eigengenes with traits
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Save correlations
write.csv(moduleTraitCor, file.path(output_dir, "module_trait_correlations.csv"))
write.csv(moduleTraitPvalue, file.path(output_dir, "module_trait_pvalues.csv"))

# Visualize module-trait relationships
pdf(file.path(output_dir, "module_trait_relationships.pdf"), width = 10, height = 12)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))

dev.off()

# ============================================================================
# 7. IDENTIFY KEY MODULES FOR EACH TRAIT
# ============================================================================

# Function to find genes in modules of interest
get_module_genes <- function(module_color, geneModulesMembership, moduleTraitCor, trait_name) {
  module <- paste0("ME", module_color)
  
  if (module %in% rownames(moduleTraitCor)) {
    module_genes <- modules %>%
      filter(module_color == !!module_color) %>%
      pull(gene_id)
    
    # Get gene significance and module membership
    gs <- as.data.frame(cor(datExpr[, module_genes], traitData[[trait_name]], use = "p"))
    colnames(gs) <- "GeneSignificance"
    gs$gene_id <- module_genes
    
    mm <- as.data.frame(cor(datExpr[, module_genes], MEs[[module]], use = "p"))
    colnames(mm) <- "ModuleMembership"
    mm$gene_id <- module_genes
    
    result <- merge(gs, mm, by = "gene_id")
    result$module <- module_color
    
    return(result)
  } else {
    return(NULL)
  }
}

# Find modules most correlated with each trait
for (trait in names(traitData)) {
  cat("\nAnalyzing trait:", trait, "\n")
  
  # Get top correlated modules
  trait_cors <- moduleTraitCor[, trait]
  top_modules <- names(sort(abs(trait_cors), decreasing = TRUE)[1:5])
  top_modules <- gsub("ME", "", top_modules)
  
  cat("Top 5 modules for", trait, ":\n")
  for (mod in top_modules) {
    cor_val <- moduleTraitCor[paste0("ME", mod), trait]
    pval <- moduleTraitPvalue[paste0("ME", mod), trait]
    cat(sprintf("  %s: r=%.3f, p=%.3e\n", mod, cor_val, pval))
  }
  
  # Save top module genes
  all_module_genes <- list()
  for (mod in top_modules) {
    module_genes <- get_module_genes(mod, MEs, moduleTraitCor, trait)
    if (!is.null(module_genes)) {
      all_module_genes[[mod]] <- module_genes
      
      # Save individual module
      write.csv(module_genes,
                file.path(output_dir, paste0("module_", mod, "_genes_", trait, ".csv")),
                row.names = FALSE)
    }
  }
}

# ============================================================================
# 8. VISUALIZE MODULE EIGENGENES
# ============================================================================
pdf(file.path(output_dir, "module_eigengenes.pdf"), width = 12, height = 10) # Height 10 is sufficient for multi-page

# 1. Prepare Data & Fix Sample ID Mismatch
ME_df <- as.data.frame(MEs)
ME_df$sample_id <- rownames(ME_df)

# FIX: R replaces hyphens with dots in column names during import. 
# We must ensure the metadata sample_ids match the WGCNA sample_ids (dots)
# Check if WGCNA IDs have dots while metadata has something else
if(any(grepl("\\.", ME_df$sample_id)) && !any(grepl("\\.", sample_info$sample_id))) {
  print("Detected mismatch: Converting metadata sample_ids to match WGCNA format (dots).")
  # Create a temporary column for merging that matches R's 'check.names' format
  sample_info$merge_id <- make.names(sample_info$sample_id)
} else {
  sample_info$merge_id <- sample_info$sample_id
}

# Merge using the corrected ID
ME_df <- merge(ME_df, sample_info, by.x = "sample_id", by.y = "merge_id")

# Debugging: Check if merge worked
if(nrow(ME_df) == 0) {
  stop("Merge failed! No common samples found. Check rownames(MEs) vs sample_info$sample_id")
}
if(is.null(ME_df$condition)) {
  stop("Column 'condition' not found after merge! Check your sample_metadata.csv")
}

# 2. Plotting with Pagination (Fixes 'Margins too large' error)
module_names <- names(MEs)

# We use a fixed layout (e.g., 6 plots per page). 
# PDF will automatically create new pages when we fill one up.
plots_per_page <- 6
par(mfrow = c(3, 2)) # 3 rows, 2 columns

for (module in module_names) {
  module_color <- gsub("ME", "", module)
  
  # Ensure the condition is a factor for proper boxplot ordering
  ME_df$condition <- as.factor(ME_df$condition)
  
  boxplot(ME_df[[module]] ~ ME_df$condition,
          main = paste("Module", module_color),
          xlab = "",
          ylab = "Eigengene Expression",
          col = module_color,
          las = 2,          # Rotate x-axis labels vertical
          outline = FALSE)  # Optional: Hide outliers for cleaner view
}

dev.off()

# ============================================================================
# 9. TARGETED MODULE ENRICHMENT ANALYSIS
# ============================================================================

# Set output directory for these specific results
module_enrich_dir <- file.path(output_dir, "Targeted_Module_Enrichment")
dir.create(module_enrich_dir, showWarnings = FALSE)

# ============================================================================
# 9.1. DEFINE ENRICHMENT FUNCTIONS
# ============================================================================

# Function for GO Enrichment
run_module_GO <- function(gene_list, module_name) {
  
  # Convert Symbols to Entrez
  gene.df <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)
  
  if (nrow(gene.df) == 0) {
    cat("  No valid Entrez IDs found for", module_name, "\n")
    return(NULL)
  }
  
  # Run GO (Biological Process)
  ego <- enrichGO(gene = gene.df$ENTREZID,
                  #universe = universe_entrez,
                  OrgDb = org.Dr.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  # Save and Plot if results exist
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), 
              file.path(module_enrich_dir, paste0("GO_BP_", module_name, ".csv")), 
              row.names = FALSE)
    
    pdf(file.path(module_enrich_dir, paste0("GO_Plot_", module_name, ".pdf")), width = 10, height = 8)
    print(dotplot(ego, showCategory = 20, title = paste("GO Enrichment:", module_name, "Module")))
    print(barplot(ego, showCategory = 20, title = paste("GO Enrichment:", module_name, "Module")))
    dev.off()
    return(ego)
  } else {
    cat("  No significant GO terms found for", module_name, "\n")
    return(NULL)
  }
}

# Function for KEGG Enrichment
run_module_KEGG <- function(gene_list, module_name) {
  
  # Convert Symbols to Entrez
  gene.df <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)
  
  if (nrow(gene.df) == 0) return(NULL)
  
  # Run KEGG
  kk <- enrichKEGG(gene = gene.df$ENTREZID,
                   #universe = universe_entrez,
                   organism = 'dre',
                   pvalueCutoff = 0.05)
  
  # Save and Plot
  if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    kk <- setReadable(kk, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")
    
    write.csv(as.data.frame(kk), 
              file.path(module_enrich_dir, paste0("KEGG_", module_name, ".csv")), 
              row.names = FALSE)
    
    pdf(file.path(module_enrich_dir, paste0("KEGG_Plot_", module_name, ".pdf")), width = 10, height = 8)
    print(dotplot(kk, showCategory = 20, title = paste("KEGG Pathway:", module_name, "Module")))
    dev.off()
    return(kk)
  } else {
    cat("  No significant KEGG pathways found for", module_name, "\n")
    return(NULL)
  }
}

# ============================================================================
# 9.2. EXTRACT GENES AND RUN ANALYSIS
# ============================================================================
#"Enrichment analysis was performed using the standard organism-wide background to ensure broad capture of functional categories."
# Define the modules of interest based on your biological question
all_modules <- gsub("ME", "", names(MEs))

cat("Starting targeted enrichment analysis...\n")

for (color in all_modules) {
  
  cat("\nProcessing Module:", color, "...\n")
  
  # A. Extract Gene IDs (Symbols) for this module color
  # Note: colnames(datExpr) contains your gene symbols
  module_genes <- colnames(datExpr)[moduleColors == color]
  
  cat("  Number of genes:", length(module_genes), "\n")
  
  # Save the raw gene list for reference
  write.csv(data.frame(gene_id = module_genes), 
            file.path(module_enrich_dir, paste0("GeneList_", color, ".csv")), 
            row.names = FALSE)
  
  # B. Run GO Enrichment
  go_res <- run_module_GO(module_genes, paste0(color, "_Module"))
  
  # C. Run KEGG Enrichment
  kegg_res <- run_module_KEGG(module_genes, paste0(color, "_Module"))
  
}

cat("\n=== Targeted Analysis Complete! ===\n")
cat("Results saved in:", module_enrich_dir, "\n")

# ============================================================================
# 10. WGCNA-GSEA ANALYSIS: Ranking by Module Membership (kME)
# ============================================================================

# -------------------------------------------------------------------------
# 1. PREPARE THE DATA
# -------------------------------------------------------------------------

# Load your FULL dataset if possible (not just the 5000 filtered genes)
# If you don't have the full 20k loaded, use 'datExpr' (it will work, just less power)
if (exists("vsd_mat")) {
  # If you have the full VST matrix (20k genes)
  gsea_datExpr <- t(vsd_mat) 
} else {
  # Fallback to the WGCNA filtered dataset
  gsea_datExpr <- datExpr 
}

# Ensure Eigengenes (MEs) are calculated
if (!exists("MEs")) {
  MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
}

# -------------------------------------------------------------------------
# 2. DEFINE GSEA FUNCTION
# -------------------------------------------------------------------------

run_wgcna_gsea <- function(module_color, expression_data, eigengenes, output_path) {
  
  cat(paste0("\nRunning GSEA for Module: ", module_color, "...\n"))
  
  # A. Calculate kME (Correlation between Genes and the Module Eigengene)
  target_ME <- eigengenes[[paste0("ME", module_color)]]
  
  # Calculate correlation for ALL genes in the dataset against this module
  kME_values <- cor(expression_data, target_ME, use = "p")
  
  # B. Prepare Ranked List
  # Create a data frame
  gene_rank_df <- data.frame(
    Symbol = rownames(kME_values),
    kME = kME_values[,1]
  )
  
  # Map Symbols to Entrez IDs (Required for GSEA)
  gene_map <- bitr(gene_rank_df$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dr.eg.db)
  
  # Merge and clean
  gene_rank_df <- merge(gene_rank_df, gene_map, by.x="Symbol", by.y="SYMBOL")
  
  # Create the named vector required by clusterProfiler
  # Format: Name = EntrezID, Value = kME
  geneList <- gene_rank_df$kME
  names(geneList) <- gene_rank_df$ENTREZID
  
  # SORT DESCENDING (Crucial step!)
  geneList <- sort(geneList, decreasing = TRUE)
  
  # C. Run GSEA (GO Biological Process)
  # Note: We use pvalueCutoff = 1.0 initially to capture everything, then filter
  gse_res <- gseGO(geneList     = geneList,
                   OrgDb        = org.Dr.eg.db,
                   ont          = "BP",
                   keyType      = "ENTREZID",
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE,
                   seed         = 123) # Set seed for reproducibility
  
  # D. Save Results
  if (!is.null(gse_res) && nrow(gse_res) > 0) {
    
    # Readable (Entrez -> Symbol)
    gse_res <- setReadable(gse_res, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")
    
    # Save CSV
    write.csv(as.data.frame(gse_res), 
              file.path(output_path, paste0("GSEA_GO_", module_color, ".csv")), 
              row.names = FALSE)
    
    # Plot 1: Dotplot
    pdf(file.path(output_path, paste0("GSEA_Dotplot_", module_color, ".pdf")), width=10, height=12)
    print(dotplot(gse_res, showCategory=15, split=".sign") + facet_grid(.~.sign) +
            ggtitle(paste("GSEA GO:", module_color, "Module")))
    dev.off()
    
    # Plot 2: Ridgeplot (Shows distribution of kME for pathways)
    pdf(file.path(output_path, paste0("GSEA_Ridgeplot_", module_color, ".pdf")), width=10, height=8)
    print(ridgeplot(gse_res, showCategory=15) + labs(x = "kME (Module Membership)"))
    dev.off()
    
    cat(paste("  Saved GSEA results for", module_color, "\n"))
    
  } else {
    cat(paste("  No significant GSEA terms found for", module_color, "\n"))
  }
}

# -------------------------------------------------------------------------
# 3. RUN THE LOOP
# -------------------------------------------------------------------------

# Set your output directory
gsea_dir <- file.path(output_dir, "GSEA_Analysis")
dir.create(gsea_dir, showWarnings = FALSE)

# Modules to analyze
all_modules <- gsub("ME", "", names(MEs))

for (mod in all_modules) {
  run_wgcna_gsea(mod, gsea_datExpr, MEs, gsea_dir)
}

cat("\n=== GSEA Analysis Complete ===\n")


# ============================================================================
# 11. COMPARECLUSTERS - KEGG
# ============================================================================

# 1. Setup Gene Clusters (Using the 4 Key Modules)
# ============================================================
# Extract genes for your 4 key modules
gene_clusters <- list(
  "Yellow\n(1116 UP)"   = colnames(datExpr)[moduleColors == "yellow"],
  "Green\n(1116 DOWN)"  = colnames(datExpr)[moduleColors == "green"],
  "Black\n(1245 UP)"    = colnames(datExpr)[moduleColors == "black"],
  "Turquoise\n(1245 DOWN)" = colnames(datExpr)[moduleColors == "turquoise"],
  "Blue" = colnames(datExpr)[moduleColors == "blue"],
  "Grey" = colnames(datExpr)[moduleColors == "grey"]
)

# Convert to Entrez IDs
gene_clusters_entrez <- lapply(gene_clusters, function(x) {
  bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dr.eg.db)$ENTREZID
})

# 2. Run CompareCluster
# ============================================================
cat("Running compareCluster...\n")
ck <- compareCluster(geneCluster = gene_clusters_entrez, 
                     fun = "enrichKEGG", 
                     organism = "dre",
                     pvalueCutoff = 0.05) # Keep this relaxed to find your terms

# Convert IDs to readable symbols
ck <- setReadable(ck, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")

# 1. CREATE A FILTERED VERSION OF THE RESULT
# We manually go into the object and keep only the top 5 rows (by p.adjust) for each Cluster.
ck_top5 <- ck  # Create a copy so you don't lose the original

# Filter the internal data frame
ck_top5@compareClusterResult <- ck@compareClusterResult %>%
  group_by(Cluster) %>%
  arrange(-Count) %>%        # Sort by Count (highest count first)
  slice_head(n = 5) %>%        # Take the top 5
  ungroup()

# 2. PLOT THE FILTERED OBJECT
pdf(file.path(output_dir, "enrichKEGG_CompareCluster_Top5_6modules.pdf"), width = 10, height = 8)

dotplot(ck_top5) + 
  ggtitle("Five Enriched Pathways per Module") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"),
        axis.text.y = element_text(size = 10))

dev.off()

