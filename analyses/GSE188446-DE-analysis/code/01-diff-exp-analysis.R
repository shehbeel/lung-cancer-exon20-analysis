# Script to differential expression analysis of HER20
# Author: Shehbeel Arif

# Load libraries
library(dplyr)
library(DESeq2)

################################################################################
## SET DIRECTORIES
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "01-diff-exp-analysis")


# Set output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

# Make output directories if they don't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Set input directory
data_dir <- file.path(root_dir, "data")

# Declare input file paths
metadata_file <- file.path(data_dir, "GSE188446_meta.txt")
data_file <- file.path(data_dir, "GSE188446_counts.txt")

#######
# Import metadata and data 
metadata <- readr::read_delim(metadata_file)
expression_df <- readr::read_delim(data_file)

#######
## PREPROCESS THE DATA
expression_df <- expression_df %>%
  tibble::column_to_rownames("#GENE")

# Check if the expression and metadata are in the same order
all.equal(colnames(expression_df), metadata$Sample)

# Make cluster a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    Treatment = factor(Treatment, levels = c("Vehicle", "TDM1", "TAK", "Combo", "Resistance"))
  )

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)


## Create DESeq2Dataset
# round all expression counts (if there are decimal values present in data)
gene_matrix <- round(filtered_expression_df)

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~Treatment
)

## Run Differential Expression Analysis
deseq_object <- DESeq(ddset)
# Extract results table
deseq_results <- results(deseq_object)
resultsNames(deseq_object)

# Use lfcShrink() function to obtain shrunken log fold change estimates based on 
# negative binomial distribution. This will add the estimates to your results table. 
# Using lfcShrink() can help decrease noise and preserve large differences between 
# groups (it requires that apeglm package be installed) (Zhu et al., Bioinformatics 2018).
# deseq_results <- lfcShrink(
#   deseq_object, # The original DESeq2 object after running DESeq()
#   coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
#   res = deseq_results # The original DESeq2 results table
# )

# same as above but with lfc shrinkage
# "Treatment_TDM1_vs_Vehicle"       "Treatment_TAK_vs_Vehicle"       
# "Treatment_Combo_vs_Vehicle"      "Treatment_Resistance_vs_Vehicle"
# generate results table for A vs ctl
Treatment_TDM1_vs_Vehicle_results <- results(deseq_object, name="Treatment_TDM1_vs_Vehicle")
# Another way:
# Treatment_TDM1_vs_Vehicle_results <- results(deseq_object, c("Treatment", "Vehicle", "TDM1"))
Treatment_TDM1_vs_Vehicle_results <- lfcShrink(deseq_object, 
                                               coef = "Treatment_TDM1_vs_Vehicle", 
                                               type = 'apeglm', 
                                               res = Treatment_TDM1_vs_Vehicle_results)
# Sort and filter DESeq2 results table and convert to dataframe
TDM1_df <- Treatment_TDM1_vs_Vehicle_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  TDM1_df,
  file.path(
    results_dir,
    "TDM1_vs_Vehicle_DEG.csv"
  )
)

Treatment_TAK_vs_Vehicle_results <- results(deseq_object, name="Treatment_TAK_vs_Vehicle")
Treatment_TAK_vs_Vehicle_results <- lfcShrink(deseq_object, 
                                               coef = "Treatment_TAK_vs_Vehicle", 
                                               type = 'apeglm', 
                                               res = Treatment_TAK_vs_Vehicle_results)
# Sort and filter DESeq2 results table and convert to dataframe
TAK_df <- Treatment_TAK_vs_Vehicle_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  TAK_df,
  file.path(
    results_dir,
    "TAK_vs_Vehicle_DEG.csv"
  )
)

Treatment_Combo_vs_Vehicle_results <- results(deseq_object, name="Treatment_Combo_vs_Vehicle")
Treatment_Combo_vs_Vehicle_results <- lfcShrink(deseq_object, 
                                               coef = "Treatment_Combo_vs_Vehicle", 
                                               type = 'apeglm', 
                                               res = Treatment_Combo_vs_Vehicle_results)
# Sort and filter DESeq2 results table and convert to dataframe
Combo_df <- Treatment_Combo_vs_Vehicle_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  Combo_df,
  file.path(
    results_dir,
    "Combo_vs_Vehicle_DEG.csv"
  )
)

Treatment_Resistance_vs_Vehicle_results <- results(deseq_object, name="Treatment_Resistance_vs_Vehicle")
Treatment_Resistance_vs_Vehicle_results <- lfcShrink(deseq_object, 
                                               coef = "Treatment_Resistance_vs_Vehicle", 
                                               type = 'apeglm', 
                                               res = Treatment_Resistance_vs_Vehicle_results)
Resistance_df <- Treatment_Resistance_vs_Vehicle_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  Resistance_df,
  file.path(
    results_dir,
    "Resistance_vs_Vehicle_DEG.csv"
  )
)




Treatment_Resistance_vs_Combo_results <- results(deseq_object, c("Treatment", "Combo", "Resistance"))
# Treatment_TDM1_vs_Vehicle_results <- lfcShrink(deseq_object, 
#                                                coef = "Treatment_Resistance_vs_Combo", 
#                                                type = 'apeglm', 
#                                                res = Treatment_Resistance_vs_Combo_results)

# Sort and filter DESeq2 results table and convert to dataframe
Resistance_vs_Combo_df <- Treatment_Resistance_vs_Combo_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  Resistance_vs_Combo_df,
  file.path(
    results_dir,
    "Resistance_vs_Combo_DEG.csv"
  )
)

Treatment_Resistance_vs_TAK_results <- results(deseq_object, c("Treatment", "TAK", "Resistance"))

# Sort and filter DESeq2 results table and convert to dataframe
Resistance_vs_TAK_df <- Treatment_Resistance_vs_TAK_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  Resistance_vs_TAK_df,
  file.path(
    results_dir,
    "Resistance_vs_TAK_DEG.csv"
  )
)

Treatment_Resistance_vs_TDM1_results <- results(deseq_object, c("Treatment", "TDM1", "Resistance"))

# Sort and filter DESeq2 results table and convert to dataframe
Resistance_vs_TDM1_df <- Treatment_Resistance_vs_TDM1_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  Resistance_vs_TDM1_df,
  file.path(
    results_dir,
    "Resistance_vs_TDM1_DEG.csv"
  )
)



################################################################################
## Create volcano plot
TDM1_vs_Vehicle_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  TDM1_df,
  lab = TDM1_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
TDM1_vs_Vehicle_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = TDM1_vs_Vehicle_volcano_plot,
  file.path(plots_dir, "TDM1_vs_Vehicle_volcano_plot.png"),
  scale=1.5
)

TAK_vs_Vehicle_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  TAK_df,
  lab = TAK_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
TAK_vs_Vehicle_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = TAK_vs_Vehicle_volcano_plot,
  file.path(plots_dir, "TAK_vs_Vehicle_volcano_plot.png"),
  scale=1.5
)

Combo_vs_Vehicle_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  Combo_df,
  lab = Combo_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
Combo_vs_Vehicle_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = Combo_vs_Vehicle_volcano_plot,
  file.path(plots_dir, "Combo_vs_Vehicle_volcano_plot.png"),
  scale=1.5
)

Resistance_vs_Vehicle_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  Resistance_df,
  lab = Resistance_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
Resistance_vs_Vehicle_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = Resistance_vs_Vehicle_volcano_plot,
  file.path(plots_dir, "Resistance_vs_Vehicle_volcano_plot.png"),
  scale=1.5
)

Resistance_vs_Combo_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  Resistance_df,
  lab = Resistance_vs_Combo_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
Resistance_vs_Combo_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = Resistance_vs_Combo_volcano_plot,
  file.path(plots_dir, "Resistance_vs_Combo_volcano_plot.png"),
  scale=1.5
)

Resistance_vs_TAK_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  Resistance_df,
  lab = Resistance_vs_TAK_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
Resistance_vs_TAK_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = Resistance_vs_TAK_volcano_plot,
  file.path(plots_dir, "Resistance_vs_TAK_volcano_plot.png"),
  scale=1.5
)

Resistance_vs_TDM1_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  Resistance_df,
  lab = Resistance_vs_TDM1_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
Resistance_vs_TDM1_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = Resistance_vs_TDM1_volcano_plot,
  file.path(plots_dir, "Resistance_vs_TDM1_volcano_plot.png"),
  scale=1.5
)
