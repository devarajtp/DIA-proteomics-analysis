# DIA-proteomics-analysis
# ============================================================================
# Paired Differential Expression Analysis - Simple Linear Model
# Filtering: Keep any protein detected at least once in 6 samples
# ============================================================================

# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(limma)
  library(ggplot2)
})

# 2. SET PATHS AND WORKING DIRECTORY
working_dir <- "fole path"
if (!dir.exists(working_dir)) dir.create(working_dir, recursive = TRUE)
setwd(working_dir)

data_path <- "mydata.xlsx"
raw_data <- read_excel(data_path)

# 3. DATA CLEANING & METADATA
all_intensity_cols <- colnames(raw_data)[grepl("^exp_1__", colnames(raw_data))]
for (col in all_intensity_cols) {
  raw_data[[col]] <- suppressWarnings(as.numeric(raw_data[[col]]))
}

# Define the 3 successful pairs (Excluding 2 and 3)
keep_pairs <- c(1, 4, 5)
sample_info <- data.frame(sample_name = all_intensity_cols, stringsAsFactors = FALSE)
sample_info$condition <- ifelse(grepl("CTRL", sample_info$sample_name), "CTRL", "Hrs_KD")
sample_info$pair_id <- as.numeric(sub(".*\\s(\\d+)$", "\\1", sample_info$sample_name))

sample_info <- sample_info[sample_info$pair_id %in% keep_pairs, ]
intensity_cols_filtered <- sample_info$sample_name

# 4. PROTEIN INFO & RAW MATRIX
proteins_info <- data.frame(
  Protein_Accession = raw_data$`Protein Accession`,
  Protein_Name = raw_data$`Protein Name`,
  Gene_ID = raw_data$`Gene ID`,
  Protein_Description = raw_data$`Protein Description`,
  stringsAsFactors = FALSE
)

expr_matrix <- as.matrix(raw_data[, intensity_cols_filtered])
rownames(expr_matrix) <- proteins_info$Protein_Accession

# 5. STEP 5: MAXIMAL FILTERING & VALID VALUE COUNTING
print("\n=== STEP 5: MAXIMAL FILTERING & COUNTING ===")

ctrl_cols <- sample_info$sample_name[sample_info$condition == "CTRL"]
kd_cols   <- sample_info$sample_name[sample_info$condition == "Hrs_KD"]

# Count non-NA values per protein
valid_count_ctrl <- rowSums(!is.na(expr_matrix[, ctrl_cols]))
valid_count_kd   <- rowSums(!is.na(expr_matrix[, kd_cols]))
total_valid      <- valid_count_ctrl + valid_count_kd

# Rule: Keep if detected in at least 1 of 6 samples
keep_idx <- total_valid >= 1

expr_filtered     <- expr_matrix[keep_idx, ]
proteins_filtered <- proteins_info[keep_idx, ]

# Store counts for the final table
counts_to_add <- data.frame(
  Protein_Accession = rownames(expr_filtered),
  Valid_Values_CTRL = valid_count_ctrl[keep_idx],
  Valid_Values_KD   = valid_count_kd[keep_idx],
  Total_Valid       = total_valid[keep_idx]
)

print(paste("Proteins retained:", nrow(expr_filtered)))

# 6. STEP 6: RAW LOG2 TRANSFORMATION
expr_log2 <- log2(expr_filtered + 1)

# 7. STEP 7: IMPUTATION (Down-shifted Normal)
set.seed(123)
expr_imputed <- expr_log2
for(i in 1:ncol(expr_imputed)) {
  na_mask <- is.na(expr_imputed[,i])
  col_mean <- mean(expr_imputed[,i], na.rm = TRUE)
  col_sd <- sd(expr_imputed[,i], na.rm = TRUE)
  expr_imputed[na_mask, i] <- rnorm(sum(na_mask), mean = col_mean - 1.8*col_sd, sd = 0.3*col_sd)
}

# ============================================================================
# STEP 6.5: QUALITY CONTROL - BOXPLOTS (Post-Log2, Pre-Imputation)
# ============================================================================
print("\n=== STEP 6.5: GENERATING BOXPLOT ===")

library(reshape2)
# Convert matrix to long format for ggplot
boxplot_data <- melt(expr_log2)
colnames(boxplot_data) <- c("Protein", "Sample", "Log2_Intensity")

# Merge with condition info for coloring
boxplot_data <- merge(boxplot_data, sample_info, by.x = "Sample", by.y = "sample_name")

pdf("03a_QC_Boxplot_Distributions.pdf", width = 8, height = 6)
p_box <- ggplot(boxplot_data, aes(x = Sample, y = Log2_Intensity, fill = condition)) +
  geom_boxplot(na.rm = TRUE, outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Log2 Intensity Distribution (Pairs 1, 4, 5)",
       subtitle = "Before Imputation",
       y = "Log2 Intensity",
       x = "") +
  scale_fill_manual(values = c("CTRL" = "grey", "Hrs_KD" = "orange"))

print(p_box)
dev.off()
# ============================================================================
# STEP 7.5: SAVE INTERMEDIATE DATA (Filtered and Imputed)
# ============================================================================
print("\n=== STEP 7.5: SAVING INTERMEDIATE DATA ===")

# 1. Save Filtered Raw Data (Values before log2 and imputation)
# We merge with protein info so you know which row is which
filtered_export <- cbind(proteins_filtered, expr_filtered)
write.csv(filtered_export, "01_Data_Filtered_Raw.csv", row.names = FALSE)

# 2. Save Imputed Data (The final log2 matrix used for Limma)
# These are the values the linear model actually "sees"
imputed_export <- cbind(proteins_filtered, as.data.frame(expr_imputed))
write.csv(imputed_export, "02_Data_Log2_Imputed.csv", row.names = FALSE)

print("Saved: 01_Data_Filtered_Raw.csv and 02_Data_Log2_Imputed.csv")

# ============================================================================
# STEP 7.6: QUALITY CONTROL - PCA PLOT
# ============================================================================
print("\n=== STEP 7.6: GENERATING PCA PLOT ===")

# Prepare data for PCA (Transpose so samples are rows)
pca_data <- t(expr_imputed)
pca_result <- prcomp(pca_data, scale. = TRUE)

# Calculate variance explained
percent_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2))

# Create data frame for plotting
pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Condition = sample_info$condition,
  Pair = as.character(sample_info$pair_id)
)

# Generate Plot
pdf("03_QC_PCA_Plot.pdf", width = 7, height = 6)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, shape = Pair)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text(aes(label = Pair), vjust = -1, show.legend = FALSE) + # Label points with Pair ID
  labs(title = "PCA Plot: Proteomics Samples (Pairs 1, 4, 5)",
       subtitle = "Based on Log2-Imputed Data",
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  scale_color_manual(values = c("CTRL" = "#333333", "Hrs_KD" = "#E69F00"))

print(p)
dev.off()

print("PCA Plot saved: 03_QC_PCA_Plot.pdf")

# 8. STEP 8: PAIRED LINEAR MODEL
pair <- factor(sample_info$pair_id)
condition <- factor(sample_info$condition, levels = c("CTRL", "Hrs_KD"))
design <- model.matrix(~ pair + condition)
colnames(design) <- c("Intercept", "Pair4", "Pair5", "HrsKD_vs_CTRL")

fit <- lmFit(expr_imputed, design)
fit2 <- eBayes(fit)

# 9. STEP 9: EXTRACT & MERGE RESULTS
results_all <- topTable(fit2, coef = "HrsKD_vs_CTRL", number = Inf, adjust.method = "BH")
results_all$Protein_Accession <- rownames(results_all)

# Merge with protein info AND the valid value counts
results_final <- results_all %>%
  merge(proteins_filtered, by = "Protein_Accession") %>%
  merge(counts_to_add, by = "Protein_Accession")

# 10. STEP 10: TOP 100 LISTS
top_100_down <- results_final %>% filter(logFC < 0) %>% arrange(P.Value) %>% head(100)
top_100_up   <- results_final %>% filter(logFC > 0) %>% arrange(P.Value) %>% head(100)

write.csv(top_100_down, "Top_100_Downregulated_with_Counts.csv", row.names = FALSE)
write.csv(top_100_up, "Top_100_Upregulated_with_Counts.csv", row.names = FALSE)
write.csv(results_final, "Full_Results_Simple_Model.csv", row.names = FALSE)

# 11. STEP 11: VOLCANO PLOT
fc_cutoff <- 0.58 # 1.5-fold
p_cutoff  <- 0.05

results_final$Significance <- "Not Significant"
results_final$Significance[results_final$logFC >= fc_cutoff & results_final$P.Value < p_cutoff] <- "Up"
results_final$Significance[results_final$logFC <= -fc_cutoff & results_final$P.Value < p_cutoff] <- "Down"

pdf("Volcano_Plot.pdf", width = 8, height = 6)
ggplot(results_final, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Down" = "blue", "Up" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  labs(title = "Volcano Plot (1.5-fold cutoff)", x = "log2 Fold Change", y = "-log10 P-value")
dev.off()

print("Analysis Complete. Files saved to working directory.")

# ============================================================================
# STEP 14: SAVE MODEL METADATA
# ============================================================================
print("\n=== STEP 14: SAVING METADATA ===")

# 1. Save the Sample Key (Mapping sample names to Condition and Rat Pair)
write.csv(sample_info, "00_Metadata_Sample_Key.csv", row.names = FALSE)

# 2. Save the Design Matrix (The actual math layout used by Limma)
write.csv(as.data.frame(design), "00_Metadata_Model_Design_Matrix.csv", row.names = TRUE)

print("Saved: 00_Metadata_Sample_Key.csv and 00_Metadata_Model_Design_Matrix.csv")
