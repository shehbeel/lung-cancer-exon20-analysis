# Gene Expression Boxplots
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia
# Source: https://rdrr.io/github/dviraran/xCell/src/R/xCell.R

## LOAD LIBRARIES
library(dplyr)
library(ggpubr)

## SET DIRECTORIES
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "GSE188446-DE-analysis")

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
mrna_data_file <- file.path(data_dir, "GSE188446_counts.txt")

#######
# Import metadata and data 
metadata <- readr::read_delim(metadata_file)
expression_df <- read.table(mrna_data_file, row.names=1)

metadata$Treatment <- as.factor(metadata$Treatment)

mrna_df <- expression_df %>%
  t() %>%
  as.data.frame() %>%
  inner_join(metadata, by=c("GENE" = "sample_id"))


################################################################################
# Make box plot with stats
# Specify statistical comparisons
my_comparisons <- list(c("Vehicle","TDM1"),
                       c("Vehicle","TAK"),
                       c("Vehicle","Combo"),
                       c("TAK","Resistance"),
                       c("TAK","Combo"))
################################################################################
## pri-miRNA -> pre-miRNA
Nppa_bp <- ggboxplot(mrna_df,
                       # Specify x values
                       x = "Treatment",
                       # Specify y values
                       y = "Nppa",
                       # Color in the box plot
                       fill = "Treatment",
                       # Specify color palette
                       palette = "jco",
                       # Add x-axis label
                       xlab="Treatment",
                       # Add y-axis label
                       ylab="Nppa Gene Expression",
                       # Add points
                       add = "jitter") +
  # Add p-value
  stat_compare_means(comparisons = my_comparisons) +
  # Add global p-value
  stat_compare_means(label.y = 16000) 

# View boxplot
print(DROSHA_bp)

# Save plot
ggsave(file.path(plots_dir, "DROSHA_cluster_boxplot.tiff"), DROSHA_bp, width = 10, height = 6)


