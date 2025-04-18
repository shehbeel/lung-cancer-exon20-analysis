# Differential expression analysis of GSE181673 dataset
# Author: Shehbeel Arif

# Install required packages if not already installed
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

# Load libraries
library(dplyr)
library(biomaRt)

################################################################################
## SET DIRECTORIES
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "GSE181673-DE-analysis")

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
metadata_file <- file.path(data_dir, "GSE181673_meta.csv")
data_file <- file.path(data_dir, "GSE181673_counts.csv")

#######
# Import metadata and data 
metadata <- readr::read_delim(metadata_file)
expression_df <- readr::read_delim(data_file)

#########################################################################################
## biomart

# Assume the gene IDs are in the first column named "GENE"
ensembl_ids <- expression_df$GENE

# Set up the connection to Ensembl using the human dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve the gene symbols corresponding to the human Ensembl IDs
# For human, we use "external_gene_name" as the gene symbol attribute.
gene_conversion <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                         filters = "ensembl_gene_id",
                         values = ensembl_ids,
                         mart = ensembl)

# Merge the conversion table with your count data
merged_data <- merge(gene_conversion, expression_df, by.x = "ensembl_gene_id", by.y = "GENE")

# Replace blank or missing external_gene_name values with the ensembl_gene_id
merged_data <- merged_data %>%
  mutate(external_gene_name = ifelse(is.na(external_gene_name) | external_gene_name == "",
                                     ensembl_gene_id,
                                     external_gene_name))

# Save results as CSV
readr::write_csv(
  merged_data,
  file.path(results_dir, "GSE181673_counts_gene_symbol.csv")
)


# ##################################################################################
# ##Alternative Using org.Hs.eg.db
# 
# # Install Bioconductor packages if necessary
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#   BiocManager::install("org.Hs.eg.db")
# }
# if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
#   BiocManager::install("AnnotationDbi")
# }
# 
# library(org.Hs.eg.db)
# library(AnnotationDbi)
# 
# # Read your CSV count data
# count_data <- read.csv("GSE181673_counts.csv", header = TRUE, stringsAsFactors = FALSE)
# 
# # Assume the gene IDs are in the first column named "GENE"
# ensembl_ids <- count_data$GENE
# 
# # Map the human Ensembl IDs to gene symbols using the org.Hs.eg.db package
# gene_symbols <- mapIds(org.Hs.eg.db,
#                        keys = ensembl_ids,
#                        column = "SYMBOL",
#                        keytype = "ENSEMBL",
#                        multiVals = "first")
# 
# # Add the gene symbols as a new column to your data frame
# count_data$GeneSymbol <- gene_symbols
# 
# # Display the first few rows of the updated data frame
# head(count_data)
# 
