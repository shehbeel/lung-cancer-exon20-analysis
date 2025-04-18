# GSVA GSE181673 - MitoCarta3.0 All Pathways GSVA
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_03_gsva.html)

## LOAD LIBRARIES
# Attach the DESeq2 library
library(DESeq2)
# Attach the `qusage` library
#library(qusage)
# Attach the `GSVA` library
library(GSVA)
# Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)
# We will need this so we can use the pipe: %>%
library(magrittr)
library(dplyr)
library(ggplot2)



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
data_file <- file.path(data_dir, "GSE181673_counts_gene_symbol.csv")
# Customs Geneset for GSEA
pathways_file <- file.path(data_dir, "Human_Genes_MitoCarta3.csv")

#######
# Import metadata and data 
metadata <- readr::read_csv(metadata_file)
expression_df <- readr::read_csv(data_file)

library(data.table)
pathways_gs <- fread(pathways_file)

#######
## PREPROCESS THE DATA
expression_df <- expression_df %>%
  dplyr::select(-c(ensembl_gene_id)) %>%
  distinct(external_gene_name, .keep_all = TRUE) %>% # Remove duplicate genes
  tibble::column_to_rownames("external_gene_name")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$Sample_ID)

# Check if the expression and metadata are in the same order
all.equal(colnames(expression_df), metadata$Sample_ID)

#############################
# Define a minimum counts cutoff and filter the data to include
expression_df <- expression_df %>%
  # Only keep rows that have total counts above the cutoff
  dplyr::filter(rowSums(.) >= 50) %>%
  # The next DESeq2 functions need the values to be converted to integers
  round()

# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = expression_df, # Our prepped data frame with counts
  colData = metadata, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)

# Perform DESeq2 normalization and transformation
# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)

# Retrieve the normalized data from the `DESeqDataSet`
vst_df <- assay(dds_norm) %>%
  as.data.frame() %>% # Make into a data frame
  tibble::rownames_to_column("gene_symbol") # Make Gene IDs into their own column

############################
# Remove duplicate rows based on a specific column (e.g., "column_name")
# mito_pathways_gs_unique <- mito_pathways_gs %>%
#   distinct(entrez_gene, .keep_all = TRUE)

# Modify the column to keep only values before ">"
pathways_gs$Minor_Pathways <- sub("\\|.*", "", pathways_gs$MitoCarta3.0_MitoPathways)


# Split the column by the ">" delimiter (trimming optional surrounding spaces)
library(stringr)
pathways_gs <- pathways_gs %>%
  mutate(MitoPathways_list = str_split(`MitoCarta3.0_MitoPathways`, "\\s*>\\s*"))






# Modify the column to keep only values before the second ">"
pathways_gs$Minor_Pathways <- sub("^([^>]*>[^>]*)>.*", "\\1", pathways_gs$MitoCarta3.0_MitoPathways)
# Modify the column to keep only values before "|"
pathways_gs$Minor_Pathways2 <- sub("\\|.*", "", pathways_gs$Minor_Pathways)
# Remove trailing spaces from the column
pathways_gs$Minor_Pathways2 <- sub("\\s+$", "", pathways_gs$Minor_Pathways2)

# Make gene sets into a list, which is required for GSVA
pathways_list <- split(
  pathways_gs$HumanGeneID, # The genes we want split into pathways
  pathways_gs$Minor_Pathways2 # The pathways made as the higher levels of the list
)

# convert entrez gene id from integer to character type as required by GSVA
#pathways_list <- lapply(pathways_list, function(genes) as.character(genes))

# Look at the gene set list
head(pathways_list, n = 2)

# GENE IDENTIFIER CONVERSION
#keytypes(org.Hs.eg.db)

# First let's create a mapped data frame we can join to the gene expression values
mapped_df <- data.frame(
  "entrez_id" = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = vst_df$gene_symbol,
    # Replace with the type of gene identifiers in your data
    keytype = "SYMBOL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
  # drop that from the data frame
  dplyr::filter(!is.na(entrez_id)) %>%
  # Make an `Ensembl` column to store the row names
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(vst_df, by = c("Ensembl" = "gene_symbol"))

# Count up how many Entrez IDs mapped to multiple Ensembl IDs
sum(duplicated(mapped_df$entrez_id))

# Handling duplicate gene identifiers
# First let's determine the gene means
gene_means <- rowMeans(mapped_df %>% dplyr::select(-Ensembl, -entrez_id))

# Let's add this as a column in our `mapped_df`.
mapped_df <- mapped_df %>%
  # Add gene_means as a column called gene_means
  dplyr::mutate(gene_means) %>%
  # Reorder the columns so `gene_means` column is upfront
  dplyr::select(Ensembl, entrez_id, gene_means, dplyr::everything())

filtered_mapped_df <- mapped_df %>%
  # Sort so that the highest mean expression values are at the top
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(entrez_id, .keep_all = TRUE)

# Check for duplicates again
sum(duplicated(filtered_mapped_df$entrez_id))


## Prepare data matrix for GSVA
filtered_mapped_matrix <- filtered_mapped_df %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Ensembl, -gene_means) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("entrez_id") %>%
  # Now we can convert our object into a matrix
  as.matrix()

#################################################################
## PERFORM GENE SET VARIATION ANALYSIS
# # R can often read in data from a URL
# hallmarks_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.1/h.all.v7.1.symbols.gmt"
# # QuSAGE is another pathway analysis method, the qusage package has a function
# # for reading GMT files and turning them into a list
# hallmarks_list <- qusage::read.gmt(hallmarks_url)

# Perform GSVA
gsva_results <- gsva(
  gsvaParam(
    exprData = filtered_mapped_matrix,
    geneSets = pathways_list,
    # Minimum gene set size
    minSize = 15, # v5 = 15, v6 = 1
    # Maximum gene set size
    maxSize = 500,
    # Kernel for estimation
    # Gaussian for our transformed data on log2-like scale
    kcdf = "Gaussian",
    # enrichment score is the difference between largest positive and negative
    maxDiff = TRUE
  ),
  # Use 4 cores for multiprocessing
  BPPARAM = BiocParallel::MulticoreParam(4)
)

# Print first 10 rows of GSVA results
head(gsva_results[, 1:10])

# Write GSVA results as csv
gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  readr::write_csv(file.path(
    results_dir,
    "MitoCarta3_allpathways_gsva_results.csv"
  ))


pheatmap::pheatmap(gsva_results, cluster_cols = FALSE)
################################################################################
# Make Complex heatmap for six clusters

# Set up column annotation from metadata
col_annot_df <- metadata %>%
  # Only select the columns of interest
  #dplyr::select(Sample_ID, Tissue, Sample_Type, Disease, Condition, Treatment) %>%
  dplyr::select(Sample_ID, Cancer, Tissue, Drug, Treatment_Time) %>%
  # Arrange by cluster and short_histology
  dplyr::arrange(Cancer, Tissue, Drug, Treatment_Time) %>%
  # Store sample
  tibble::column_to_rownames("Sample_ID")


# Create the ComplexHeatmap column annotation object
col_annot <- ComplexHeatmap::HeatmapAnnotation(
  # Supply Treatment labels
  Cancer = col_annot_df$Cancer,
  Tissue = col_annot_df$Tissue,
  Drug = col_annot_df$Drug,
  Treatment_Time = col_annot_df$Treatment_Time,
  # Pick colors for each Treatment
  col = list(Cancer = c("NSCLC" = "#114B5F",
                        "Breast Carcinoma" = "#1A936F",
                        "Esophagus Carcinoma" = "#88D498"
                        ),
             
            Tissue = c("AU565 cancer cell line" = "#071E22", 
                        "CTG-2534 PDX model" = "#1D7874", 
                        "NCI-H1781 cancer cell line" = "#679289", 
                        "NCI-H2170-YVMA CRISPR/Cas9 engineered cell line" = "#F4C095", 
                        "OE19 cancer cell line" = "#EE2E31"
                        ),
            
            Drug = c("DMSO" = "#C9CBA3",
                     "Natrosol" = "#FFE1A8",
                     "BI-1622" = "#E26D5C",
                     "BI-4142" = "#723D46",
                     "Poziotinib" = "#472D30"
                        ),
            Treatment_Time = c(
                      "0h" = "#FAF3DD",
                      "0.5h" = "#C8D5B9",
                      "2h" = "#8FC0A9",
                      "6h" = "#68B0AB",
                      "16h" = "#4A7C59",
                      "24h" = "#364958"
                        )
  )
)


###########################


# Create a color function based on standardized scale
color_func <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("navy", "white", "firebrick3")
)

# Check if column annotation df and expression data are in the same order
all.equal(colnames(gsva_results), rownames(col_annot_df))
# Re-order data
gsva_results <- gsva_results %>%
  as.data.frame() %>%
  dplyr::select(rownames(col_annot_df)) %>%
  as.matrix()

# custom_order <- c("Innate Immune", "Innate Immune: Canonical", "Innate Immune: Non-Canonical", "Innate Immune: Inflammation",                            
#                   "Adaptive Immune", "Adaptive Immune: Surface Marker/Receptor Signaling", "mtDNA/dsRNA Release", "mtDNA/dsRNA Innate Immune",                              
#                   "mtDNA/dsRNA Innate Immune: Activated by mtDNA/dsRNA", "ISR", "ISR: Target-Gene", "ISR: Target-Gene: Cytokine/Chemokine",                   
#                   "UPR MT/ER", "UPR MT/ER: Target-Gene", "UPR MT/ER: Target-Gene: UPR ER", "PANoptosis",                                             
#                   "OXPHOS", "OXPHOS: Complex I", "OXPHOS: Complex III", "OXPHOS: Complex IV", "OXPHOS: Complex V", "OXPHOS: MT-Ribosome",                                    
#                   "Mitochondrial Metabolism", "Mitochondrial Metabolism: Glycolysis",                   
#                   "Mitochondrial Metabolism: TCA Cycle", "Mitochondrial Metabolism: Fatty Acid Oxidation",         
#                   "Mitochondrial Metabolism: Nucleotide Synthesis", "Mitochondrial Metabolism: Nucleotide Synthesis: Nuclear",
#                   "Antioxidant Defenses", "Cytosolic Protein Import", "HIF", "HIF: Selected HIF target", "mTOR", "mTOR: mTOR Structural", "mTOR: mTOR Structural: mTORC1", "Peroxisome")
# 
# # Reorder rows based on custom_order
# gsva_results <- gsva_results[custom_order, ]









# Color annotation bar on left/right showing main pathway groups
# add row gaps between main groups
# column gaps between cell line groups

## Separate heatmaps for each cell line

# Generate col_split using the Comparison column
column_split <- factor(col_annot_df$Cancer, levels = unique(col_annot_df$Cancer))


gsva_hm <- ComplexHeatmap::Heatmap(gsva_results, 
                                   name = "GSVA score",
                                   show_row_names = TRUE,
                                   row_names_side = "left",
                                   show_column_names = FALSE,
                                   cluster_rows = FALSE, 
                                   cluster_columns = FALSE,
                                   show_row_dend = FALSE,
                                   column_split = column_split,
                                   top_annotation = col_annot,
                                   row_names_gp = gpar(fontsize = 10),  # Reduce font size for row names
                                   col = color_func#,
                                   #width = unit(20, "cm"),
                                   #height = unit(20, "cm")
)


library(ComplexHeatmap)
# fig1_heatmap <- h1 %v% h2
gsva_heatmap <- ggplotify::as.ggplot(grid.grabExpr(
  draw(gsva_hm, padding = unit(c(5, 50, 5, 5), "mm")) # Add extra space for row names
))

# Save the Heatmap plot
ggsave(
  plot = gsva_heatmap,
  file.path(plots_dir, "MitoCarta3.0_allpathways_gsva_heatmap.png"),
  width = 20,
  height = 6
)
