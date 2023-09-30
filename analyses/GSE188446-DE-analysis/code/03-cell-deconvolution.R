# Cell Deconvolution using Bulk RNA-Seq Data (xCell)
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia
# Source: https://rdrr.io/github/dviraran/xCell/src/R/xCell.R

## LOAD LIBRARIES
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)

## SET DIRECTORIES
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "GSE188446-DE-analysis")
# Source this script, which contains a wrapper function that to calculate correlations
source(file.path(root_dir, "utils", "xCell_function.R"))

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
mrna_data_file <- file.path(data_dir, "GSE188446_TPM.txt")
xCell.data.file <- file.path(data_dir, "xCell.data.rda")

#######
# Load expression and xCell data
metadata <- readr::read_tsv(metadata_file)
# Expression data needs to have gene symbols as row names and columns as samples
mrna_expression_df <- readr::read_tsv(mrna_data_file) #%>% 
#tibble::column_to_rownames(var="Geneid")
load(xCell.data.file)
#######
# Map Ensembl IDs to Gene Symbols
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

# First let's create a mapped data frame we can join to the gene expression values
mapped_df <- data.frame(
  "GENE" = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = mrna_expression_df$Geneid,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "SYMBOL",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Gene Symbol,
  # drop that from the data frame
  dplyr::filter(!is.na(GENE)) %>%
  # Make an `Ensembl` column to store the row names
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(mrna_expression_df, by = c("Ensembl" = "Geneid"))

# Count up how many Entrez IDs mapped to multiple Ensembl IDs
sum(duplicated(mapped_df$GENE))

# Handling duplicate gene identifiers
# First let's determine the gene means
gene_means <- rowMeans(mapped_df %>% dplyr::select(-Ensembl, -GENE))

# Let's add this as a column in our `mapped_df`.
mapped_df <- mapped_df %>%
  # Add gene_means as a column called gene_means
  dplyr::mutate(gene_means) %>%
  # Reorder the columns so `gene_means` column is upfront
  dplyr::select(Ensembl, GENE, gene_means, dplyr::everything())

filtered_mapped_df <- mapped_df %>%
  # Sort so that the highest mean expression values are at the top
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(GENE, .keep_all = TRUE)

# Check for duplicates again
sum(duplicated(filtered_mapped_df$GENE))

# GSVA can't the Ensembl IDs so we should drop this column as well as the means
filtered_mapped_df <- filtered_mapped_df %>%
  dplyr::select(-Ensembl, -gene_means) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("GENE") #%>%
# Now we can convert our object into a matrix
#as.matrix()

#######
# Run xCell
xCell_results <- xCellAnalysis(filtered_mapped_df) 
xCell_results <- xCell_results %>%
  as.data.frame()

# Save output xCell Results
readr::write_csv(xCell_results, file = file.path(results_dir, ".PBTA_miRNA_xCell_Deconvolution.csv"))

# Plot xCell results Heatmap
pheatmap::pheatmap(scale(xCell_results))



# Create a color function based on standardized scale
color_func <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#67a9cf", "#f7f7f7", "#ef8a62")
  #c("blue", "white", "red")
)

# Set up column annotation from metadata
col_annot_df <- metadata %>%
  # Only select the treatment and sample ID columns
  dplyr::select(Sample_ID, Treatment, Tissue) %>%
  # Add on the eigengene expression by joining with sample IDs
  # dplyr::inner_join(module_eigengene, by = "refinebio_accession_code") %>%
  # Arrange by cluster and short_histology
  dplyr::arrange(Treatment, Tissue) %>%
  # Store sample
  tibble::column_to_rownames("Sample_ID")

col_annot <- ComplexHeatmap::HeatmapAnnotation(
  # Supply Column labels
  Treatment = col_annot_df$Treatment,
  Tissue = col_annot_df$Tissue#,
  # Pick colors for each miRNA cluster in cluster
  #col = list(cluster = c("1" = "#0073C2FF", "2" = "#EFC000FF", "3" = "#868686FF", "4" = "#CD534CFF"))
)

xcell_heatmap <- ComplexHeatmap::Heatmap(as.matrix(scale(xCell_results)), 
                                         name = "z-score",
                                         show_row_names = TRUE,
                                         row_names_side = "right",
                                         show_column_names = FALSE,
                                         cluster_rows = TRUE, 
                                         cluster_columns = FALSE,
                                         top_annotation = col_annot,
                                         col = color_func#,
                                         #width = unit(20, "cm"),
                                         #height = unit(20, "cm")
)

xcell_heatmap <- ggplotify::as.ggplot(grid.grabExpr(draw(xcell_heatmap)))

# Save volcano plot
ggsave(
  plot = xcell_heatmap,
  file.path(plots_dir, "xCell_Heatmap.png"),
  width = 16,
  height = 14
)

