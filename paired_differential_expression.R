# ============================================================================
# Paired Differential Expression Analysis - Simple Linear Model
# Filtering: Keep any protein detected at least once in all samples
# ============================================================================

# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(limma)
  library(ggplot2)
  library(reshape2)
})

# ============================================================================
# USER CONFIGURATION - Edit these parameters before running the script
# ============================================================================

# Path to your input Excel file (use forward slashes or double backslashes)
data_path <- "data.xlsx"

# Name of the working/output directory (will be created if it does not exist)
working_dir <- "R_Analysis_Output"

# Prefix used to identify intensity columns in your data
intensity_col_prefix <- "exp_1__"

# Condition labels as they appear in the intensity column names
condition_1_label <- "CTRL"     # e.g., control group
condition_2_label <- "KD"       # e.g., knockdown / treatment group

# Pairs to include in the analysis (set to NULL to include all detected pairs)
# Example: keep_pairs <- c(1, 4, 5)  to include only pairs 1, 4 and 5
# Set to NULL to include all pairs automatically
keep_pairs <- NULL

# Minimum number of samples a protein must be detected in (non-NA) to be retained
min_valid_samples <- 1

# Fold change cutoff for the volcano plot (in log2 scale; 0.58 ≈ 1.5-fold)
fc_cutoff <- 0.58

# P-value cutoff for the volcano plot
p_cutoff <- 0.05

# ============================================================================
# END OF USER CONFIGURATION
# ============================================================================

# 2. SET WORKING DIRECTORY AND RESOLVE DATA PATH
# Resolve data_path to an absolute path before changing the working directory
if (!file.exists(data_path)) stop("Input file not found: ", data_path)
data_path_abs <- normalizePath(data_path)

if (!dir.exists(working_dir)) dir.create(working_dir, recursive = TRUE)
setwd(working_dir)

# 3. LOAD DATA
raw_data <- read_excel(data_path_abs)

# 4. DATA CLEANING & METADATA
all_intensity_cols <- colnames(raw_data)[grepl(paste0("^", intensity_col_prefix), colnames(raw_data))]
for (col in all_intensity_cols) {
  raw_data[[col]] <- suppressWarnings(as.numeric(raw_data[[col]]))
}

sample_info <- data.frame(sample_name = all_intensity_cols, stringsAsFactors = FALSE)
sample_info$condition <- ifelse(
  grepl(condition_1_label, sample_info$sample_name),
  condition_1_label,
  condition_2_label
)
sample_info$pair_id <- as.numeric(sub(".*\\s(\\d+)$", "\\1", sample_info$sample_name))

# Apply pair filtering if specified
if (!is.null(keep_pairs)) {
  sample_info <- sample_info[sample_info$pair_id %in% keep_pairs, ]
}
intensity_cols_filtered <- sample_info$sample_name

# 5. PROTEIN INFO & RAW MATRIX
proteins_info <- data.frame(
  Protein_Accession = raw_data$`Protein Accession`,
  Protein_Name = raw_data$`Protein Name`,
  Gene_ID = raw_data$`Gene ID`,
  Protein_Description = raw_data$`Protein Description`,
  stringsAsFactors = FALSE
)

expr_matrix <- as.matrix(raw_data[, intensity_cols_filtered])
rownames(expr_matrix) <- proteins_info$Protein_Accession

# 6. FILTERING & VALID VALUE COUNTING
print("\n=== FILTERING & COUNTING ===")

condition1_cols <- sample_info$sample_name[sample_info$condition == condition_1_label]
condition2_cols <- sample_info$sample_name[sample_info$condition == condition_2_label]

valid_count_cond1 <- rowSums(!is.na(expr_matrix[, condition1_cols]))
valid_count_cond2 <- rowSums(!is.na(expr_matrix[, condition2_cols]))
total_valid       <- valid_count_cond1 + valid_count_cond2

keep_idx <- total_valid >= min_valid_samples

expr_filtered     <- expr_matrix[keep_idx, ]
proteins_filtered <- proteins_info[keep_idx, ]

counts_to_add <- data.frame(
  Protein_Accession       = rownames(expr_filtered),
  Valid_Values_Condition1 = valid_count_cond1[keep_idx],
  Valid_Values_Condition2 = valid_count_cond2[keep_idx],
  Total_Valid             = total_valid[keep_idx]
)

print(paste("Proteins retained:", nrow(expr_filtered)))

# 7. LOG2 TRANSFORMATION
expr_log2 <- log2(expr_filtered + 1)

# 8. IMPUTATION (Down-shifted Normal)
set.seed(123)
expr_imputed <- expr_log2
for (i in 1:ncol(expr_imputed)) {
  na_mask <- is.na(expr_imputed[, i])
  col_mean <- mean(expr_imputed[, i], na.rm = TRUE)
  col_sd   <- sd(expr_imputed[, i], na.rm = TRUE)
  expr_imputed[na_mask, i] <- rnorm(sum(na_mask), mean = col_mean - 1.8 * col_sd, sd = 0.3 * col_sd)
}

# ============================================================================
# QUALITY CONTROL - BOXPLOTS (Post-Log2, Pre-Imputation)
# ============================================================================
print("\n=== GENERATING BOXPLOT ===")

boxplot_data <- melt(expr_log2)
colnames(boxplot_data) <- c("Protein", "Sample", "Log2_Intensity")
boxplot_data <- merge(boxplot_data, sample_info, by.x = "Sample", by.y = "sample_name")

pdf("03a_QC_Boxplot_Distributions.pdf", width = 8, height = 6)
p_box <- ggplot(boxplot_data, aes(x = Sample, y = Log2_Intensity, fill = condition)) +
  geom_boxplot(na.rm = TRUE, outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Log2 Intensity Distribution",
       subtitle = "Before Imputation",
       y = "Log2 Intensity",
       x = "") +
  scale_fill_manual(values = setNames(c("grey", "orange"), c(condition_1_label, condition_2_label)))
print(p_box)
dev.off()

# ============================================================================
# SAVE INTERMEDIATE DATA (Filtered and Imputed)
# ============================================================================
print("\n=== SAVING INTERMEDIATE DATA ===")

filtered_export <- cbind(proteins_filtered, expr_filtered)
write.csv(filtered_export, "01_Data_Filtered_Raw.csv", row.names = FALSE)

imputed_export <- cbind(proteins_filtered, as.data.frame(expr_imputed))
write.csv(imputed_export, "02_Data_Log2_Imputed.csv", row.names = FALSE)

print("Saved: 01_Data_Filtered_Raw.csv and 02_Data_Log2_Imputed.csv")

# ============================================================================
# QUALITY CONTROL - PCA PLOT
# ============================================================================
print("\n=== GENERATING PCA PLOT ===")

pca_data   <- t(expr_imputed)
pca_result <- prcomp(pca_data, scale. = TRUE)
percent_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2))

pca_df <- data.frame(
  PC1       = pca_result$x[, 1],
  PC2       = pca_result$x[, 2],
  Condition = sample_info$condition,
  Pair      = as.character(sample_info$pair_id)
)

pdf("03_QC_PCA_Plot.pdf", width = 7, height = 6)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, shape = Pair)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text(aes(label = Pair), vjust = -1, show.legend = FALSE) +
  labs(title = "PCA Plot: Proteomics Samples",
       subtitle = "Based on Log2-Imputed Data",
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  scale_color_manual(values = setNames(c("#333333", "#E69F00"), c(condition_1_label, condition_2_label)))
print(p)
dev.off()

print("PCA Plot saved: 03_QC_PCA_Plot.pdf")

# 9. PAIRED LINEAR MODEL
pair      <- factor(sample_info$pair_id)
condition <- factor(sample_info$condition, levels = c(condition_1_label, condition_2_label))
design    <- model.matrix(~ pair + condition)

# Rename the condition coefficient for readability
coef_name <- paste0(condition_2_label, "_vs_", condition_1_label)
colnames(design)[ncol(design)] <- coef_name

fit  <- lmFit(expr_imputed, design)
fit2 <- eBayes(fit)

# 10. EXTRACT & MERGE RESULTS
results_all <- topTable(fit2, coef = coef_name, number = Inf, adjust.method = "BH")
results_all$Protein_Accession <- rownames(results_all)

results_final <- results_all %>%
  merge(proteins_filtered, by = "Protein_Accession") %>%
  merge(counts_to_add, by = "Protein_Accession")

# 11. TOP 100 LISTS
top_100_down <- results_final %>% filter(logFC < 0) %>% arrange(P.Value) %>% head(100)
top_100_up   <- results_final %>% filter(logFC > 0) %>% arrange(P.Value) %>% head(100)

write.csv(top_100_down, "Top_100_Downregulated_with_Counts.csv", row.names = FALSE)
write.csv(top_100_up,   "Top_100_Upregulated_with_Counts.csv",   row.names = FALSE)
write.csv(results_final, "Full_Results_Simple_Model.csv",         row.names = FALSE)

# 12. VOLCANO PLOT
results_final$Significance <- "Not Significant"
results_final$Significance[results_final$logFC >= fc_cutoff  & results_final$P.Value < p_cutoff] <- "Up"
results_final$Significance[results_final$logFC <= -fc_cutoff & results_final$P.Value < p_cutoff] <- "Down"

pdf("Volcano_Plot.pdf", width = 8, height = 6)
ggplot(results_final, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Down" = "blue", "Up" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  labs(title = paste0("Volcano Plot (", round(2^fc_cutoff, 2), "-fold cutoff)"),
       x = "log2 Fold Change",
       y = "-log10 P-value")
dev.off()

# 13. SAVE MODEL METADATA
print("\n=== SAVING METADATA ===")

write.csv(sample_info, "00_Metadata_Sample_Key.csv", row.names = FALSE)
write.csv(as.data.frame(design), "00_Metadata_Model_Design_Matrix.csv", row.names = TRUE)

print("Saved: 00_Metadata_Sample_Key.csv and 00_Metadata_Model_Design_Matrix.csv")
print("Analysis Complete. Files saved to working directory.")
