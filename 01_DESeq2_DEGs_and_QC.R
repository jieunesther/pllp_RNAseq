# ============================================================================
# Comprehensive RNA-seq Analysis Pipeline
# Experiment: pllp mutants (EK, 1116, 1245) x Diet (HP, LP)
# ============================================================================

# Load required libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(VennDiagram)
library(ggrepel)
library(edgeR)
library(limma)
library(WGCNA)
library(GO.db)
library(org.Dr.eg.db)
library(clusterProfiler)
library(enrichplot)
library(grid)

# Create output directory
dir.create("../analysis_results", showWarnings = FALSE)
output_dir <- "../analysis_results"

# ============================================================================
# 1. SAMPLE METADATA AND DATA IMPORT
# ============================================================================

# Create sample metadata
sample_info <- data.frame(
  sample_id = c(
    "EK_HP_1_S1", "EK_HP_3_S2", "EK_HP_4_S3",
    "EK_LP_5_S4", "EK_LP_7_S5", "EK_LP_8_S6",
    "pd1116_HP_9_S7", "pd1116_HP_10_S8", "pd1116_HP_11_S9",
    "pd1116_LP_13_S10", "pd1116_LP_14_S11", "pd1116_LP_15_S12",
    "pd1245_HP_17_S13", "pd1245_HP_18_S14", "pd1245_HP_19_S15",
    "pd1245_LP_21_S16", "pd1245_LP_23_S17", "pd1245_LP_25_S18"
  ),
  genotype = rep(c("EK", "EK", "pd1116", "pd1116", "pd1245", "pd1245"), each = 3),
  diet = rep(c("HP", "LP"), each = 3, times = 3),
  replicate = rep(1:3, times = 6)
)

# Add filename column
sample_info$filename <- paste0(sample_info$sample_id, ".count")

# Convert to factors with appropriate reference levels
sample_info$genotype <- factor(sample_info$genotype, levels = c("EK", "pd1116", "pd1245"))
sample_info$diet <- factor(sample_info$diet, levels = c("HP", "LP"))

# Create combined condition column for pairwise comparisons
sample_info$condition <- factor(paste(sample_info$genotype, sample_info$diet, sep = "_"))

print("Sample metadata:")
print(sample_info)
write.csv(sample_info, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)

# ============================================================================
# 2. IMPORT COUNT DATA
# ============================================================================

# Import HTSeq count files
count_data <- NULL
for(i in 1:nrow(sample_info)) {
  counts <- read.table(sample_info$filename[i], 
                       header = FALSE, 
                       stringsAsFactors = FALSE,
                       col.names = c("gene_id", sample_info$sample_id[i]))
  
  if(is.null(count_data)) {
    count_data <- counts
  } else {
    count_data <- merge(count_data, counts, by = "gene_id")
  }
}

# Remove HTSeq summary rows
count_data <- count_data %>%
  filter(!grepl("^__", gene_id))

# Set gene IDs as rownames
rownames(count_data) <- count_data$gene_id
count_data <- count_data[, -1]

# Ensure columns are in same order as sample_info
count_data <- count_data[, sample_info$sample_id]

print(paste("Total genes:", nrow(count_data)))
print(paste("Total samples:", ncol(count_data)))

# Save raw counts
write.csv(count_data, file.path(output_dir, "raw_counts.csv"), quote = FALSE)

# ============================================================================
# 3. QUALITY CONTROL AND FILTERING
# ============================================================================

# Calculate total counts per sample
sample_info$total_counts <- colSums(count_data)

# Plot library sizes
pdf(file.path(output_dir, "QC_library_sizes.pdf"), width = 10, height = 6)
ggplot(sample_info, aes(x = sample_id, y = total_counts/1e6, fill = condition)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Total Reads (millions)", title = "Library Sizes") +
  scale_fill_brewer(palette = "Set3")
dev.off()

# Filter low count genes (keep genes with at least 10 counts in at least 3 samples)
keep <- rowSums(count_data >= 10) >= 3
count_data_filtered <- count_data[keep, ]
print(paste("Genes after filtering:", nrow(count_data_filtered)))

# Function to convert gene IDs to symbols (if available)
get_gene_symbols <- function(gene_ids) {
  tryCatch({
    # Try to get gene symbols from org.Dr.eg.db
    library(org.Dr.eg.db)
    
    # Convert ENSEMBL to SYMBOL
    symbols <- mapIds(org.Dr.eg.db,
                      keys = gene_ids,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
    
    # If symbol is NA, keep the original ID
    symbols <- ifelse(is.na(symbols), gene_ids, symbols)
    
    # Make unique (in case of duplicate symbols)
    symbols <- make.unique(as.character(symbols))
    
    return(symbols)
    
  }, error = function(e) {
    # If conversion fails, return original IDs
    cat("Gene symbol conversion failed, using gene IDs\n")
    return(gene_ids)
  })
}

gene_symbol <- get_gene_symbols(rownames(count_data_filtered))
rownames(count_data_filtered) <- gene_symbol

# ============================================================================
# 4. DESeq2 ANALYSIS
# ============================================================================

# Create DESeq2 object with interaction model
dds <- DESeqDataSetFromMatrix(
  countData = count_data_filtered,
  colData = sample_info,
  design = ~ genotype + diet + genotype:diet
)

# Run DESeq2
dds <- DESeq(dds)

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file.path(output_dir, "normalized_counts.csv"), quote = FALSE)

# Variance stabilizing transformation for visualization
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)
write.csv(vsd_mat, file.path(output_dir, "vst_transformed_counts.csv"), quote = FALSE)

# rlog transformation (alternative to VST)
rld <- rlog(dds, blind = FALSE)
rld_mat <- assay(rld)

# ============================================================================
# 5. PRINCIPAL COMPONENT ANALYSIS (PCA)
# ============================================================================

# PCA using VST data
pdf(file.path(output_dir, "PCA_all_samples.pdf"), width = 10, height = 7)

# PCA colored by genotype and diet
pcaData <- plotPCA(vsd, intgroup = c("genotype", "diet"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

my_colors <- c("EK" = "#E8A87C", "pd1116" = "#8CD493", "pd1245" = "#789BCA")

p1 <- ggplot(pcaData, aes(PC1, PC2, color=genotype, shape=diet)) +
  geom_point(size=5) + # Increased point size slightly for better visibility
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  #ggtitle("PCA Analysis: Genotype and Diet") +
  scale_color_manual(values = my_colors) +
  
  # --- CUSTOM FONT SIZES START HERE ---
  theme(
    # 1. Axis Titles (PC1: X% variance, etc.)
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    
    # 2. Axis Text (The numbers on the axis)
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    
    # 3. Legend Title (Genotype, Diet)
    legend.title = element_text(size = 14, face = "bold"),
    
    # 4. Legend Text (EK, mut1116, HP, LP)
    legend.text = element_text(size = 12),
    
    # 5. Plot Title
    #plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
print(p1)

dev.off()

# 3D PCA plot data
pca_result <- prcomp(t(vsd_mat))
pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  PC3 = pca_result$x[,3],
  sample_info
)
write.csv(pca_df, file.path(output_dir, "PCA_coordinates.csv"), row.names = FALSE)

# ============================================================================
# 6. DIFFERENTIAL EXPRESSION ANALYSIS - ALL COMPARISONS
# ============================================================================

# Function to extract and save results
extract_results <- function(dds, contrast, name, alpha = 0.05, lfc_threshold = 1) {
  res <- results(dds, contrast = contrast, alpha = alpha)
  res_ordered <- res[order(res$padj), ]
  
  # Add gene symbols if available
  res_df <- as.data.frame(res_ordered)
  res_df$gene_id <- rownames(res_df)
  
  # Write all results
  write.csv(res_df, 
            file.path(output_dir, paste0("DESeq2_results_", name, ".csv")),
            row.names = FALSE)
  
  # Extract significant genes
  sig_genes <- res_df %>%
    filter(padj < alpha & abs(log2FoldChange) > lfc_threshold)
  
  write.csv(sig_genes,
            file.path(output_dir, paste0("DESeq2_significant_", name, ".csv")),
            row.names = FALSE)
  
  # Summary
  cat("\n", name, ":\n")
  cat("Total DEGs (padj < ", alpha, ", |LFC| > ", lfc_threshold, "): ", nrow(sig_genes), "\n")
  cat("Upregulated: ", sum(sig_genes$log2FoldChange > 0), "\n")
  cat("Downregulated: ", sum(sig_genes$log2FoldChange < 0), "\n")
  
  return(list(all = res_df, sig = sig_genes))
}

# Main effect comparisons
cat("\n=== MAIN EFFECTS ===\n")

# ============================================================================
# 7. CONDITION-SPECIFIC COMPARISONS
# ============================================================================

cat("\n=== CONDITION-SPECIFIC COMPARISONS ===\n")

# Recreate DESeq2 object with condition as design (no interaction term)
dds_condition <- DESeqDataSetFromMatrix(
  countData = count_data_filtered,
  colData = sample_info,
  design = ~ condition
)

dds_condition <- DESeq(dds_condition)

# All pairwise comparisons
comparisons <- list(
  # HP conditions
  list(c("condition", "pd1116_HP", "EK_HP"), "HP_1116_vs_EK"),
  list(c("condition", "pd1245_HP", "EK_HP"), "HP_1245_vs_EK"),
  list(c("condition", "pd1245_HP", "pd1116_HP"), "HP_1245_vs_1116"),
  
  # LP conditions
  list(c("condition", "pd1116_LP", "EK_LP"), "LP_1116_vs_EK"),
  list(c("condition", "pd1245_LP", "EK_LP"), "LP_1245_vs_EK"),
  list(c("condition", "pd1245_LP", "pd1116_LP"), "LP_1245_vs_1116"),
  
  # Diet effects within genotypes
  list(c("condition", "EK_HP", "EK_LP"), "EK_HP_vs_LP"),
  list(c("condition", "pd1116_HP", "pd1116_LP"), "1116_HP_vs_LP"),
  list(c("condition", "pd1245_HP", "pd1245_LP"), "1245_HP_vs_LP")
)

# Store all results
all_comparisons <- list()
for(comp in comparisons) {
  all_comparisons[[comp[[2]]]] <- extract_results(dds_condition, comp[[1]], comp[[2]])
}

# ============================================================================
# 8. VENN DIAGRAMS - DEG OVERLAPS
# ============================================================================

pdf(file.path(output_dir, "venn_diagrams.pdf"), width = 10, height = 10)

# PLOT 1: HP DIET (Genotype Effect)
# ============================================================================

# Define lists first for clarity
list_HP <- list(
  "HP_1116" = all_comparisons$HP_1116_vs_EK$sig$gene_id,
  "HP_1245" = all_comparisons$HP_1245_vs_EK$sig$gene_id
)

venn.plot1 <- venn.diagram(
  x = list_HP,
  filename = NULL,
  
  # --- 1. VISUAL SETUP ---
  main = "DEGs: Genotype Effects in HP Diet",
  fill = c("#8CD493", "#789BCA"),  
  alpha = 0.5,
  
  # --- 2. CRITICAL FIXES (Size & Order) ---
  scaled = FALSE,      # KEEPS CIRCLES SAME SIZE (prevents re-ordering based on size)
  euler.d = FALSE,     # Ensure standard Venn behavior
  
  # --- 3. LABELS (Removing them) ---
  category.names = c("", ""), # Removes the labels so you can add your own
  
  # --- 4. NUMBERS & PERCENTAGES ---
  # "raw" = count, "percent" = %, this adds the (xx%) below the number automatically
  print.mode = c("raw", "percent"), 
  
  # --- 5. FONT CONTROLS ---
  main.cex = 3,        # Title size
  cex = 3.5,           # Number size (Increased from 2.5 to 3.5)
  main.fontface = "bold",
  fontface = "bold",
  
  # Margins
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot1)


# PLOT 2: LP DIET (Genotype Effect)
# ============================================================================

list_LP <- list(
  "LP_1116" = all_comparisons$LP_1116_vs_EK$sig$gene_id,
  "LP_1245" = all_comparisons$LP_1245_vs_EK$sig$gene_id
)

venn.plot2 <- venn.diagram(
  x = list_LP,
  filename = NULL,
  
  # --- 1. VISUAL SETUP ---
  main = "DEGs: Genotype Effects in LP Diet",
  fill = c("#8CD493", "#789BCA"),
  alpha = 0.5,
  
  # --- 2. CRITICAL FIXES ---
  scaled = FALSE,      # Keeps circles constant size
  euler.d = FALSE,
  rotation.degree = 180, # <--- THE FIX: Rotates 180 degrees to put 1116 on Left, 1245 on Right
  
  # --- 3. REMOVE LABELS ---
  category.names = c("", ""),
  
  # --- 4. NUMBERS & PERCENTAGES ---
  print.mode = c("raw", "percent"),
  
  # --- 5. FONT CONTROLS ---
  main.cex = 3,
  cex = 3.5,           # Large numbers inside the circle
  main.fontface = "bold",
  fontface = "bold",
  
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot2)


# PLOT 3: DIET EFFECTS (Matched Format)
# ============================================================================

venn.plot3 <- venn.diagram(
  x = list(
    "EK"     = all_comparisons$EK_HP_vs_LP$sig$gene_id,
    "pd1116" = all_comparisons$`1116_HP_vs_LP`$sig$gene_id,
    "pd1245" = all_comparisons$`1245_HP_vs_LP`$sig$gene_id
  ),
  filename = NULL,
  
  # --- 1. VISUAL SETUP ---
  main = "DEGs: Diet Effects (LP vs HP)",
  fill = c("#E8A87C", "#8CD493", "#789BCA"), # Orange (EK), Green (1116), Blue (1245)
  alpha = 0.5,
  
  # --- 2. CRITICAL FIXES ---
  scaled = FALSE,       # Keep circles same size
  euler.d = FALSE,      # Force standard 3-circle Venn
  
  # --- 3. LABELS & NUMBERS ---
  category.names = c("", "", ""),   # Remove labels
  print.mode = c("raw", "percent"), # Show Count + Percentage
  
  # --- 4. FONT CONTROLS ---
  main.cex = 3,         # Title Size
  cex = 3,            # Number Size (Matched to others)
  main.fontface = "bold",
  fontface = "bold",
  
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot3)

dev.off()

# ============================================================================
# 9. OVERLAPPING GENES
# ============================================================================

# Save overlap analysis
overlaps <- list(
  HP_genotype_overlap = intersect(
    all_comparisons$HP_1116_vs_EK$sig$gene_id,
    all_comparisons$HP_1245_vs_EK$sig$gene_id
  ),
  LP_genotype_overlap = intersect(
    all_comparisons$LP_1116_vs_EK$sig$gene_id,
    all_comparisons$LP_1245_vs_EK$sig$gene_id
  ),
  diet_all_genotypes = Reduce(intersect, list(
    all_comparisons$EK_HP_vs_LP$sig$gene_id,
    all_comparisons$`1116_HP_vs_LP`$sig$gene_id,
    all_comparisons$`1245_HP_vs_LP`$sig$gene_id
  )),
  "1245_vs_1116_both_diets" = intersect(
    all_comparisons$HP_1245_vs_1116$sig$gene_id,
    all_comparisons$LP_1245_vs_1116$sig$gene_id
  )
)

# Genes unique to 1116 (potential compensation genes)
genes_1116_unique_HP <- setdiff(
  all_comparisons$HP_1116_vs_EK$sig$gene_id,
  all_comparisons$HP_1245_vs_EK$sig$gene_id
)

genes_1116_unique_HP_table <- all_comparisons$HP_1116_vs_EK$sig %>%
  filter(gene_id %in% genes_1116_unique_HP)

genes_1116_unique_LP <- setdiff(
  all_comparisons$LP_1116_vs_EK$sig$gene_id,
  all_comparisons$LP_1245_vs_EK$sig$gene_id
)

genes_1116_unique_LP_table <- all_comparisons$LP_1116_vs_EK$sig %>%
  filter(gene_id %in% genes_1116_unique_LP)

# Genes unique to 1245 
genes_1245_unique_HP <- setdiff(
  all_comparisons$HP_1245_vs_EK$sig$gene_id,
  all_comparisons$HP_1116_vs_EK$sig$gene_id
)

genes_1245_unique_HP_table <- all_comparisons$HP_1245_vs_EK$sig %>%
  filter(gene_id %in% genes_1245_unique_HP)

genes_1245_unique_LP <- setdiff(
  all_comparisons$LP_1245_vs_EK$sig$gene_id,
  all_comparisons$LP_1116_vs_EK$sig$gene_id
)

genes_1245_unique_LP_table <- all_comparisons$LP_1245_vs_EK$sig %>%
  filter(gene_id %in% genes_1245_unique_LP)

# Save unique genes
write.csv(
  genes_1116_unique_HP_table,
  file.path(output_dir, "genes_unique_to_1116_HP_table.csv"),
  row.names = FALSE
)

write.csv(
  genes_1116_unique_LP_table,
  file.path(output_dir, "genes_unique_to_1116_LP_table.csv"),
  row.names = FALSE
)

write.csv(
  genes_1245_unique_HP_table,
  file.path(output_dir, "genes_unique_to_1245_HP_table.csv"),
  row.names = FALSE
)

write.csv(
  genes_1245_unique_LP_table,
  file.path(output_dir, "genes_unique_to_1245_LP_table.csv"),
  row.names = FALSE
)

# ============================================================================
# 10. DEG SUMMARY PLOT & UPSET PLOT
# ============================================================================

library(ggplot2)
library(tidyverse)

# 1. Clean and Organize Data (Same as before)
plot_data <- summary_table %>%
  mutate(
    Category = case_when(
      grepl("HP vs LP", Comparison) ~ "Diet Effects",
      TRUE ~ "Genotype Effects"
    ),
    Down_Plot = -Downregulated # Negative values for left bars
  )

# Set logical order
plot_data$Comparison <- factor(plot_data$Comparison, levels = rev(plot_data$Comparison))


# 2. Generate the Plot
deg_summary <- ggplot(plot_data) +
  # --- BARS ---
  # We map 'fill' to a specific name inside aes() to force a legend to appear
  geom_col(aes(x = Comparison, y = Down_Plot, fill = "Downregulated"), width = 0.7) +
  geom_col(aes(x = Comparison, y = Upregulated, fill = "Upregulated"), width = 0.7) +
  
  # --- TEXT LABELS ---
  geom_text(aes(x = Comparison, y = Down_Plot, label = Downregulated), 
            hjust = 1.2, color = "cornflowerblue", fontface = "bold", size = 5) + # Increased size
  geom_text(aes(x = Comparison, y = Upregulated, label = Upregulated), 
            hjust = -0.2, color = "firebrick3", fontface = "bold", size = 5) +    # Increased size
  
  # --- SCALES & AXES ---
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  
  # Custom limits: -800 (Left) to +500 (Right)
  scale_y_continuous(labels = abs, limits = c(-800, 500)) + 
  
  # Manual colors for the legend
  scale_fill_manual(name = "Direction", 
                    values = c("Downregulated" = "cornflowerblue", 
                               "Upregulated" = "firebrick3")) +
  
  coord_flip() +
  
  # --- VISUAL STYLING ---
  # Increased base_size to 16 for bigger general text
  theme_minimal(base_size = 16) + 
  
  labs(
    title = "Differential Expression Summary",
    subtitle = "Count of Significant Genes (FDR < 0.05, |LFC| > 1)",
    x = "", 
    y = "Number of DEGs"
  ) +
  
  theme(
    panel.grid.major.y = element_blank(),
    
    # Axis Text
    axis.text.y = element_text(face = "bold", color = "black", size = 14),
    axis.text.x = element_text(size = 12),
    
    # Facet Labels (The Grey Boxes)
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "grey90", color = NA),
    
    # Title Centering
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 20)), # Add margin below subtitle
    
    # Legend Styling
    legend.position = "top", # Put legend at the top to save side space
    legend.title = element_blank(), # Remove the word "Direction"
    legend.text = element_text(size = 14, face = "bold")
  )

pdf(file.path(output_dir, "DEG_summary_plot.pdf"), width = 10, height = 10)
print(deg_summary)
dev.off()

library(UpSetR)
# Create a list of your gene sets
# (I am using the names from your bar plot for consistency)
upset_input_list <- list(
  "HP: pd1116 vs EK"  = all_comparisons$HP_1116_vs_EK$sig$gene_id,
  "HP: pd1245 vs EK"  = all_comparisons$HP_1245_vs_EK$sig$gene_id,
  "LP: pd1116 vs EK"  = all_comparisons$LP_1116_vs_EK$sig$gene_id,
  "LP: pd1245 vs EK"  = all_comparisons$LP_1245_vs_EK$sig$gene_id
)

upset_input_list2 <- list(
  "EK: HP vs LP"      = all_comparisons$EK_HP_vs_LP$sig$gene_id,
  "pd1116: HP vs LP"  = all_comparisons$`1116_HP_vs_LP`$sig$gene_id,
  "pd1245: HP vs LP"  = all_comparisons$`1245_HP_vs_LP`$sig$gene_id
)

upset_input_list3 <- list(
  "HP: pd1116 vs EK"  = all_comparisons$HP_1116_vs_EK$sig$gene_id,
  "HP: pd1245 vs EK"  = all_comparisons$HP_1245_vs_EK$sig$gene_id
)

upset_input_list4 <- list(
  "LP: pd1116 vs EK"  = all_comparisons$LP_1116_vs_EK$sig$gene_id,
  "LP: pd1245 vs EK"  = all_comparisons$LP_1245_vs_EK$sig$gene_id
)

pdf(file.path(output_dir, "upset_plot.pdf"), width = 10, height = 10)
upset1 <- upset(fromList(upset_input_list), 
                # 1. SHOW ALL SETS (Important! Default is only 5)
                nsets = 9, 
                
                # 2. Sort intersections by size (largest overlaps first)
                order.by = "freq", 
                
                # 3. Size Styling
                text.scale = c(2, 2, 2, 2, 2, 2), # Scale up all text elements
                point.size = 3.5,                 # Size of the dots in the matrix
                line.size = 1.5,                  # Thickness of the connecting lines
                
                # 4. Color Styling
                mainbar.y.label = "Intersection Size (Gene Count)",
                sets.x.label = "Total Set Size",
                main.bar.color = "black",
                sets.bar.color = "firebrick3",    # Color of the horizontal bars on the left
                matrix.color = "black"
)
upset2 <- upset(fromList(upset_input_list2), 
                # 1. SHOW ALL SETS (Important! Default is only 5)
                nsets = 9, 
                
                # 2. Sort intersections by size (largest overlaps first)
                order.by = "freq", 
                
                # 3. Size Styling
                text.scale = c(2, 2, 2, 2, 2, 2), # Scale up all text elements
                point.size = 3.5,                 # Size of the dots in the matrix
                line.size = 1.5,                  # Thickness of the connecting lines
                
                # 4. Color Styling
                mainbar.y.label = "Intersection Size (Gene Count)",
                sets.x.label = "Total Set Size",
                main.bar.color = "black",
                sets.bar.color = "firebrick3",    # Color of the horizontal bars on the left
                matrix.color = "black"
)
upset3 <- upset(fromList(upset_input_list3), 
                # 1. SHOW ALL SETS (Important! Default is only 5)
                nsets = 9, 
                
                # 2. Sort intersections by size (largest overlaps first)
                order.by = "freq", 
                
                # 3. Size Styling
                text.scale = c(2, 2, 2, 2, 2, 2), # Scale up all text elements
                point.size = 3.5,                 # Size of the dots in the matrix
                line.size = 1.5,                  # Thickness of the connecting lines
                
                # 4. Color Styling
                mainbar.y.label = "Intersection Size (Gene Count)",
                sets.x.label = "Total Set Size",
                main.bar.color = "black",
                sets.bar.color = "firebrick3",    # Color of the horizontal bars on the left
                matrix.color = "black"
)
upset4 <- upset(fromList(upset_input_list4), 
                # 1. SHOW ALL SETS (Important! Default is only 5)
                nsets = 9, 
                
                # 2. Sort intersections by size (largest overlaps first)
                order.by = "freq", 
                
                # 3. Size Styling
                text.scale = c(2, 2, 2, 2, 2, 2), # Scale up all text elements
                point.size = 3.5,                 # Size of the dots in the matrix
                line.size = 1.5,                  # Thickness of the connecting lines
                
                # 4. Color Styling
                mainbar.y.label = "Intersection Size (Gene Count)",
                sets.x.label = "Total Set Size",
                main.bar.color = "black",
                sets.bar.color = "firebrick3",    # Color of the horizontal bars on the left
                matrix.color = "black"
)
print(upset1)
print(upset2)
print(upset3)
print(upset4)
dev.off()

# ============================================================================
# 11. SUMMARY STATISTICS
# ============================================================================

# Create summary table
summary_table <- data.frame(
  Comparison = c(
    "HP: pd1116 vs EK", "HP: pd1245 vs EK", "HP: pd1245 vs pd1116",
    "LP: pd1116 vs EK", "LP: pd1245 vs EK", "LP: pd1245 vs pd1116",
    "EK: HP vs LP", "pd1116: HP vs LP", "pd1245: HP vs LP"
  ),
  Total_DEGs = c(
    nrow(all_comparisons$HP_1116_vs_EK$sig),
    nrow(all_comparisons$HP_1245_vs_EK$sig),
    nrow(all_comparisons$HP_1245_vs_1116$sig),
    nrow(all_comparisons$LP_1116_vs_EK$sig),
    nrow(all_comparisons$LP_1245_vs_EK$sig),
    nrow(all_comparisons$LP_1245_vs_1116$sig),
    nrow(all_comparisons$EK_HP_vs_LP$sig),
    nrow(all_comparisons$`1116_HP_vs_LP`$sig),
    nrow(all_comparisons$`1245_HP_vs_LP`$sig)
  ),
  Upregulated = c(
    sum(all_comparisons$HP_1116_vs_EK$sig$log2FoldChange > 1),
    sum(all_comparisons$HP_1245_vs_EK$sig$log2FoldChange > 1),
    sum(all_comparisons$HP_1245_vs_1116$sig$log2FoldChange > 1),
    sum(all_comparisons$LP_1116_vs_EK$sig$log2FoldChange > 1),
    sum(all_comparisons$LP_1245_vs_EK$sig$log2FoldChange > 1),
    sum(all_comparisons$LP_1245_vs_1116$sig$log2FoldChange > 1),
    sum(all_comparisons$EK_HP_vs_LP$sig$log2FoldChange > 1),
    sum(all_comparisons$`1116_HP_vs_LP`$sig$log2FoldChange > 1),
    sum(all_comparisons$`1245_HP_vs_LP`$sig$log2FoldChange > 1)
  ),
  Downregulated = c(
    sum(all_comparisons$HP_1116_vs_EK$sig$log2FoldChange < -1),
    sum(all_comparisons$HP_1245_vs_EK$sig$log2FoldChange < -1),
    sum(all_comparisons$HP_1245_vs_1116$sig$log2FoldChange < -1),
    sum(all_comparisons$LP_1116_vs_EK$sig$log2FoldChange < -1),
    sum(all_comparisons$LP_1245_vs_EK$sig$log2FoldChange < -1),
    sum(all_comparisons$LP_1245_vs_1116$sig$log2FoldChange < -1),
    sum(all_comparisons$EK_HP_vs_LP$sig$log2FoldChange < -1),
    sum(all_comparisons$`1116_HP_vs_LP`$sig$log2FoldChange < -1),
    sum(all_comparisons$`1245_HP_vs_LP`$sig$log2FoldChange < -1)
  )
)

write.csv(summary_table, file.path(output_dir, "DEG_summary_table.csv"), row.names = FALSE)
print("\n=== SUMMARY TABLE ===")
print(summary_table)

cat("\n=== Analysis Complete! ===\n")
cat("Results saved in:", output_dir, "\n")


