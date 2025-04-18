# GSVA GSE181673 - Custom Mito Pathways GSVA
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
#pathways_file <- file.path(data_dir, "custom_Mitochondria_Pathways-EntrezID.gmt")
pathways_file <- file.path(data_dir, "24_11.27_MT-SpaceBio_Pathways_Hm_Recent.gmt")


#######
# Import metadata and data 
metadata <- readr::read_csv(metadata_file)
expression_df <- readr::read_csv(data_file)
pathways_gs <- clusterProfiler::read.gmt(pathways_file)

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

# # Modify the column to keep only values before ">"
# pathways_gs$Minor_Pathways <- sub("\\|.*", "", pathways_gs$MitoCarta3.0_MitoPathways)
# 
# 
# # Split the column by the ">" delimiter (trimming optional surrounding spaces)
# library(stringr)
# pathways_gs <- pathways_gs %>%
#   mutate(MitoPathways_list = str_split(`MitoCarta3.0_MitoPathways`, "\\s*>\\s*"))
# 
# 
# 
# 
# 
# 
# # Modify the column to keep only values before the second ">"
# pathways_gs$Minor_Pathways <- sub("^([^>]*>[^>]*)>.*", "\\1", pathways_gs$MitoCarta3.0_MitoPathways)
# # Modify the column to keep only values before "|"
# pathways_gs$Minor_Pathways2 <- sub("\\|.*", "", pathways_gs$Minor_Pathways)
# # Remove trailing spaces from the column
# pathways_gs$Minor_Pathways2 <- sub("\\s+$", "", pathways_gs$Minor_Pathways2)

# Make gene sets into a list, which is required for GSVA
pathways_list <- split(
  pathways_gs$gene, # The genes we want split into pathways
  pathways_gs$term # The pathways made as the higher levels of the list
)

# convert entrez gene id from integer to character type as required by GSVA
#pathways_list <- lapply(pathways_list, function(genes) as.character(genes))

# Look at the gene set list
head(pathways_list, n = 2)

# GENE IDENTIFIER CONVERSION
#keytypes(org.Hs.eg.db)

# # First let's create a mapped data frame we can join to the gene expression values
# mapped_df <- data.frame(
#   "entrez_id" = mapIds(
#     # Replace with annotation package for the organism relevant to your data
#     org.Hs.eg.db,
#     keys = vst_df$gene_symbol,
#     # Replace with the type of gene identifiers in your data
#     keytype = "SYMBOL",
#     # Replace with the type of gene identifiers you would like to map to
#     column = "ENTREZID",
#     # This will keep only the first mapped value for each Ensembl ID
#     multiVals = "first"
#   )
# ) %>%
#   # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
#   # drop that from the data frame
#   dplyr::filter(!is.na(entrez_id)) %>%
#   # Make an `Ensembl` column to store the row names
#   tibble::rownames_to_column("Ensembl") %>%
#   # Now let's join the rest of the expression data
#   dplyr::inner_join(vst_df, by = c("Ensembl" = "gene_symbol"))
# 
# # Count up how many Entrez IDs mapped to multiple Ensembl IDs
# sum(duplicated(mapped_df$entrez_id))

# Handling duplicate gene identifiers
# First let's determine the gene means
gene_means <- rowMeans(vst_df %>% dplyr::select(-gene_symbol))

# Let's add this as a column in our `mapped_df`.
filtered_df <- vst_df %>%
  # Add gene_means as a column called gene_means
  dplyr::mutate(gene_means) %>%
  # Reorder the columns so `gene_means` column is upfront
  dplyr::select(gene_symbol, gene_means, dplyr::everything())

filtered_df <- filtered_df %>%
  # Sort so that the highest mean expression values are at the top
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(gene_symbol, .keep_all = TRUE)

# Check for duplicates again
sum(duplicated(filtered_df$gene_symbol))


## Prepare data matrix for GSVA
filtered_matrix <- filtered_df %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-gene_means) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("gene_symbol") %>%
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
    exprData = filtered_matrix,
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
    "CustomMitoSpacePathways_gsva_results.csv"
  ))


################################################################################
# Make Complex heatmap 

# Set up column annotation from metadata
col_annot_df <- metadata %>%
  # Only select the columns of interest
  #dplyr::select(Sample_ID, Tissue, Sample_Type, Disease, Condition, Treatment) %>%
  dplyr::select(Sample_ID, Cancer, Tissue, Drug, Treatment_Time) %>%
  # Arrange by cluster and short_histology
  dplyr::arrange(Cancer, Tissue, Drug, Treatment_Time) %>%
  dplyr::arrange(factor(Treatment_Time, levels = c("0h", "0.05h", "2h", "6h", "16h", "24h"))) %>%
  dplyr::arrange(factor(Drug, levels = c("DMSO", "Natrosol", "BI-1622", "BI-4142", "Poziotinib"))) %>%
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

# ROW ANNOTATION
row_annot_file <- file.path(data_dir, "24_11.27_MT-SpaceBio_Pathways_Hm_Recent_meta.csv")
row_annot_df <- readr::read_csv(row_annot_file)
row_annot_df <- row_annot_df %>%
  dplyr::filter(Pathway %in% rownames(gsva_results)) %>%
  tibble::column_to_rownames(var = "Pathway")

# tmp_df <- gsva_results %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var="Pathway") %>%
#   dplyr::left_join(row_annot_df) %>%
#   dplyr::select("Pathway", "Pathway_Hierarchy")

row_annot <- ComplexHeatmap::rowAnnotation(
  # Supply Comparison labels
  Pathway_Hierarchy = row_annot_df$Pathway_Hierarchy,
  # Pick colors for each column
  col = list(
    Pathway_Hierarchy = c(
      "Senescence"= "#a44e3f",                # rusty red (inflammaging / stress)
      "Innate Immune" = "#c24604",             # same as Immune anchor
      "Adaptive Immune" = "#d95d18",           # lighter burnt orange
      "mtDNA/dsRNA Innate Immune" = "#e1792a", # orange-brown
      "ISR" = "#705555",                       # stress gray/brown
      "UPR MT/ER" = "#b2786a",                 # pinkish brown (ER stress)
      "PANoptosis" = "#95342e",                # dark red (cell death)
      "OXPHOS" = "#5c755b",                    # olive green anchor
      "Metabolism" = "#9ba33b",                # green-yellow anchor
      "Antioxidant Defenses" = "#758d3d",      # earthy green
      "Cytosolic Protein Import" = "#6e7b4e",  # dull olive
      "Global Translation" = "#b5be4b",        # yellow-green
      "Cyto-Ribosome" = "#a6ac3e",             # similar to Metabolism
      "HIF" = "#4c5b52",                       # hypoxia dark green-gray
      "mTOR" = "#697b3c",                      # regulatory green
      "Mitophagy" = "#44553f",                 # dark moss green
      "RAAS" = "#914d1c"                       # brown-orange (vascular)
    )
  )
)

## Separate heatmaps for each cell line

# Generate col_split using the Comparison column
column_split <- factor(col_annot_df$Cancer, levels = unique(col_annot_df$Cancer))
row_split <- factor(row_annot_df$Pathway_Hierarchy, levels = unique(row_annot_df$Pathway_Hierarchy))

###########################

# Check if column annotation df and expression data are in the same order
all.equal(colnames(gsva_results), rownames(col_annot_df))
# Re-order data
gsva_results <- gsva_results %>%
  as.data.frame() %>%
  dplyr::select(rownames(col_annot_df)) %>%
  as.matrix()

# Create a color function based on standardized scale
color_func <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("navy", "white", "firebrick3")
)


library(ComplexHeatmap)
gsva_hm <- ComplexHeatmap::Heatmap(gsva_results, 
                                   name = "GSVA score",
                                   show_row_names = TRUE,
                                   row_names_side = "right",
                                   show_column_names = FALSE,
                                   cluster_rows = FALSE, 
                                   cluster_columns = FALSE,
                                   show_row_dend = FALSE,
                                   column_split = column_split,
                                   row_split = row_split,
                                   top_annotation = col_annot,
                                   left_annotation = row_annot,
                                   row_names_gp = gpar(fontsize = 10),  # Reduce font size for row names
                                   col = color_func#,
                                   #width = unit(20, "cm"),
                                   #height = unit(20, "cm")
)



# fig1_heatmap <- h1 %v% h2
gsva_heatmap <- ggplotify::as.ggplot(grid.grabExpr(
  draw(gsva_hm, padding = unit(c(5, 40, 5, 5), "mm")) # Add extra space for row names
))

# Save the Heatmap plot
ggsave(
  plot = gsva_heatmap,
  file.path(plots_dir, "CustomMitoSpacePathways_gsva_heatmap.png"),
  width = 21,
  height = 14
)


################################################################################
## ESOPHAGUS

# Set up column annotation from metadata
col_annot_df <- metadata %>%
  dplyr::filter(Cancer == "Esophagus Carcinoma") %>%
  # Only select the columns of interest
  dplyr::select(Sample_ID, Tissue, Drug, Treatment_Time) %>%
  # Arrange by cluster and short_histology
  dplyr::arrange(Drug, Treatment_Time) %>%
  dplyr::arrange(factor(Treatment_Time, levels = c("6h", "24h"))) %>%
  dplyr::arrange(factor(Drug, levels = c("DMSO", "BI-4142", "Poziotinib"))) %>%
  # Store sample
  tibble::column_to_rownames("Sample_ID")


# Create the ComplexHeatmap column annotation object
col_annot <- ComplexHeatmap::HeatmapAnnotation(
  # Supply Treatment labels
  #Cancer = col_annot_df$Cancer,
  #Tissue = col_annot_df$Tissue,
  Drug = col_annot_df$Drug,
  Treatment_Time = col_annot_df$Treatment_Time,
  # Pick colors for each Treatment
  col = list(
  # Cancer = c("NSCLC" = "#114B5F",
  #                       "Breast Carcinoma" = "#1A936F",
  #                       "Esophagus Carcinoma" = "#88D498"
  # ),
  
  # Tissue = c("AU565 cancer cell line" = "#071E22", 
  #            "CTG-2534 PDX model" = "#1D7874", 
  #            "NCI-H1781 cancer cell line" = "#679289", 
  #            "NCI-H2170-YVMA CRISPR/Cas9 engineered cell line" = "#F4C095", 
  #            "OE19 cancer cell line" = "#EE2E31"
  # ),
  
  Drug = c("DMSO" = "#C9CBA3",
           #"Natrosol" = "#FFE1A8",
           #"BI-1622" = "#E26D5C",
           "BI-4142" = "#723D46",
           "Poziotinib" = "#472D30"
  ),
  Treatment_Time = c(
    #"0h" = "#FAF3DD",
    #"0.5h" = "#C8D5B9",
    #"2h" = "#8FC0A9",
    "6h" = "#68B0AB",
    #"16h" = "#4A7C59",
    "24h" = "#364958"
  )
  )
)

# ROW ANNOTATION
row_annot_file <- file.path(data_dir, "24_11.27_MT-SpaceBio_Pathways_Hm_Recent_meta_Arif_updated.csv")
row_annot_df <- readr::read_csv(row_annot_file)
row_annot_df <- row_annot_df %>%
  dplyr::filter(Pathway %in% rownames(gsva_results)) %>%
  tibble::column_to_rownames(var = "Pathway")

row_annot <- ComplexHeatmap::rowAnnotation(
  # Supply Comparison labels
  Pathway_Hierarchy = row_annot_df$Pathway_Hierarchy,
  # Pick colors for each column
  col = list(
    Pathway_Hierarchy = c(
      "Senescence"= "#a44e3f",                # rusty red (inflammaging / stress)
      "Innate Immune" = "#c24604",             # same as Immune anchor
      "Adaptive Immune" = "#d95d18",           # lighter burnt orange
      "mtDNA/dsRNA Innate Immune" = "#e1792a", # orange-brown
      "ISR" = "#705555",                       # stress gray/brown
      "UPR MT/ER" = "#b2786a",                 # pinkish brown (ER stress)
      "PANoptosis" = "#95342e",                # dark red (cell death)
      "OXPHOS" = "#5c755b",                    # olive green anchor
      "Metabolism" = "#9ba33b",                # green-yellow anchor
      "Antioxidant Defenses" = "#758d3d",      # earthy green
      "Cytosolic Protein Import" = "#6e7b4e",  # dull olive
      "Global Translation" = "#b5be4b",        # yellow-green
      "Cyto-Ribosome" = "#a6ac3e",             # similar to Metabolism
      "HIF" = "#4c5b52",                       # hypoxia dark green-gray
      "mTOR" = "#697b3c",                      # regulatory green
      "Mitophagy" = "#44553f",                 # dark moss green
      "RAAS" = "#914d1c"                       # brown-orange (vascular)
    )
  )
)

## Separate heatmaps for each cell line

# Generate col_split using the Comparison column
column_split <- factor(col_annot_df$Drug, levels = unique(col_annot_df$Drug))
row_split <- factor(row_annot_df$Pathway_Hierarchy, levels = unique(row_annot_df$Pathway_Hierarchy))

###########################

# Check if column annotation df and expression data are in the same order
all.equal(colnames(ec_gsva_results), rownames(col_annot_df))
# Re-order data
ec_gsva_results <- gsva_results %>%
  as.data.frame() %>%
  dplyr::select(rownames(col_annot_df)) %>%
  as.matrix()

# Create a color function based on standardized scale
color_func <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("navy", "white", "firebrick3")
)


library(ComplexHeatmap)
ec_gsva_hm <- ComplexHeatmap::Heatmap(ec_gsva_results, 
                                   name = "GSVA score",
                                   show_row_names = TRUE,
                                   row_names_side = "right",
                                   show_column_names = FALSE,
                                   cluster_rows = FALSE, 
                                   cluster_columns = FALSE,
                                   show_row_dend = FALSE,
                                   column_split = column_split,
                                   row_split = row_split,
                                   top_annotation = col_annot,
                                   left_annotation = row_annot,
                                   row_names_gp = gpar(fontsize = 10),  # Reduce font size for row names
                                   col = color_func#,
                                   #width = unit(20, "cm"),
                                   #height = unit(20, "cm")
)

ec_gsva_heatmap <- ggplotify::as.ggplot(grid.grabExpr(draw(ec_gsva_hm, 
                                                               heatmap_legend_side = "bottom", 
                                                               annotation_legend_side = "bottom", 
                                                               merge_legend = FALSE, 
                                                               padding = unit(c(5, 5, 5, 25), "mm")) # Add extra space for row names
))

# Save the Heatmap plot
ggsave(
  plot = ec_gsva_heatmap,
  file.path(plots_dir, "EC_CustomMitoSpacePathways_gsva_heatmap.png"),
  width = 8,
  height = 16
)

################################################################################
## BREAST
# Set up column annotation from metadata
col_annot_df <- metadata %>%
  dplyr::filter(Cancer == "Breast Carcinoma") %>%
  # Only select the columns of interest
  dplyr::select(Sample_ID, Tissue, Drug, Treatment_Time) %>%
  # Arrange by cluster and short_histology
  dplyr::arrange(Drug, Treatment_Time) %>%
  dplyr::arrange(factor(Treatment_Time, levels = c("6h", "24h"))) %>%
  dplyr::arrange(factor(Drug, levels = c("DMSO", "BI-4142", "Poziotinib"))) %>%
  # Store sample
  tibble::column_to_rownames("Sample_ID")


# Create the ComplexHeatmap column annotation object
col_annot <- ComplexHeatmap::HeatmapAnnotation(
  # Supply Treatment labels
  #Cancer = col_annot_df$Cancer,
  #Tissue = col_annot_df$Tissue,
  Drug = col_annot_df$Drug,
  Treatment_Time = col_annot_df$Treatment_Time,
  # Pick colors for each Treatment
  col = list(
    # Cancer = c("NSCLC" = "#114B5F",
    #                       "Breast Carcinoma" = "#1A936F",
    #                       "Esophagus Carcinoma" = "#88D498"
    # ),
    
    # Tissue = c("AU565 cancer cell line" = "#071E22", 
    #            "CTG-2534 PDX model" = "#1D7874", 
    #            "NCI-H1781 cancer cell line" = "#679289", 
    #            "NCI-H2170-YVMA CRISPR/Cas9 engineered cell line" = "#F4C095", 
    #            "OE19 cancer cell line" = "#EE2E31"
    # ),
    
    Drug = c("DMSO" = "#C9CBA3",
             #"Natrosol" = "#FFE1A8",
             #"BI-1622" = "#E26D5C",
             "BI-4142" = "#723D46",
             "Poziotinib" = "#472D30"
    ),
    Treatment_Time = c(
      #"0h" = "#FAF3DD",
      #"0.5h" = "#C8D5B9",
      #"2h" = "#8FC0A9",
      "6h" = "#68B0AB",
      #"16h" = "#4A7C59",
      "24h" = "#364958"
    )
  )
)

# ROW ANNOTATION
row_annot_file <- file.path(data_dir, "24_11.27_MT-SpaceBio_Pathways_Hm_Recent_meta_Arif_updated.csv")
row_annot_df <- readr::read_csv(row_annot_file)
row_annot_df <- row_annot_df %>%
  dplyr::filter(Pathway %in% rownames(gsva_results)) %>%
  tibble::column_to_rownames(var = "Pathway")

row_annot <- ComplexHeatmap::rowAnnotation(
  # Supply Comparison labels
  Pathway_Hierarchy = row_annot_df$Pathway_Hierarchy,
  # Pick colors for each column
  col = list(
    Pathway_Hierarchy = c(
      "Senescence"= "#a44e3f",                # rusty red (inflammaging / stress)
      "Innate Immune" = "#c24604",             # same as Immune anchor
      "Adaptive Immune" = "#d95d18",           # lighter burnt orange
      "mtDNA/dsRNA Innate Immune" = "#e1792a", # orange-brown
      "ISR" = "#705555",                       # stress gray/brown
      "UPR MT/ER" = "#b2786a",                 # pinkish brown (ER stress)
      "PANoptosis" = "#95342e",                # dark red (cell death)
      "OXPHOS" = "#5c755b",                    # olive green anchor
      "Metabolism" = "#9ba33b",                # green-yellow anchor
      "Antioxidant Defenses" = "#758d3d",      # earthy green
      "Cytosolic Protein Import" = "#6e7b4e",  # dull olive
      "Global Translation" = "#b5be4b",        # yellow-green
      "Cyto-Ribosome" = "#a6ac3e",             # similar to Metabolism
      "HIF" = "#4c5b52",                       # hypoxia dark green-gray
      "mTOR" = "#697b3c",                      # regulatory green
      "Mitophagy" = "#44553f",                 # dark moss green
      "RAAS" = "#914d1c"                       # brown-orange (vascular)
    )
  )
)

## Separate heatmaps for each cell line

# Generate col_split using the Comparison column
column_split <- factor(col_annot_df$Drug, levels = unique(col_annot_df$Drug))
row_split <- factor(row_annot_df$Pathway_Hierarchy, levels = unique(row_annot_df$Pathway_Hierarchy))

###########################

# Check if column annotation df and expression data are in the same order
all.equal(colnames(bc_gsva_results), rownames(col_annot_df))
# Re-order data
bc_gsva_results <- gsva_results %>%
  as.data.frame() %>%
  dplyr::select(rownames(col_annot_df)) %>%
  as.matrix()

# Create a color function based on standardized scale
color_func <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("navy", "white", "firebrick3")
)


library(ComplexHeatmap)
bc_gsva_hm <- ComplexHeatmap::Heatmap(bc_gsva_results, 
                                      name = "GSVA score",
                                      show_row_names = TRUE,
                                      row_names_side = "right",
                                      show_column_names = FALSE,
                                      cluster_rows = FALSE, 
                                      cluster_columns = FALSE,
                                      show_row_dend = FALSE,
                                      column_split = column_split,
                                      row_split = row_split,
                                      top_annotation = col_annot,
                                      left_annotation = row_annot,
                                      row_names_gp = gpar(fontsize = 10),  # Reduce font size for row names
                                      col = color_func#,
                                      #width = unit(20, "cm"),
                                      #height = unit(20, "cm")
)

bc_gsva_heatmap <- ggplotify::as.ggplot(grid.grabExpr(draw(bc_gsva_hm, 
                                                           heatmap_legend_side = "bottom", 
                                                           annotation_legend_side = "bottom", 
                                                           merge_legend = FALSE, 
                                                           padding = unit(c(5, 5, 5, 25), "mm")) # Add extra space for row names
))

# Save the Heatmap plot
ggsave(
  plot = bc_gsva_heatmap,
  file.path(plots_dir, "BC_CustomMitoSpacePathways_gsva_heatmap.png"),
  width = 8,
  height = 16
)

################################################################################
## NSLC CTG-2534 PDX model
# Set up column annotation from metadata
col_annot_df <- metadata %>%
  dplyr::filter(Cancer == "NSCLC") %>%
  dplyr::filter(Tissue == "CTG-2534 PDX model") %>%
  # Only select the columns of interest
  dplyr::select(Sample_ID, Tissue, Drug, Treatment_Time) %>%
  # Arrange by cluster and short_histology
  dplyr::arrange(Drug, Treatment_Time) %>%
  dplyr::arrange(factor(Treatment_Time, levels = c("6h", "24h"))) %>%
  dplyr::arrange(factor(Drug, levels = c("Natrosol", "BI-4142"))) %>%
  # Store sample
  tibble::column_to_rownames("Sample_ID")


# Create the ComplexHeatmap column annotation object
col_annot <- ComplexHeatmap::HeatmapAnnotation(
  # Supply Treatment labels
  #Cancer = col_annot_df$Cancer,
  #Tissue = col_annot_df$Tissue,
  Drug = col_annot_df$Drug,
  Treatment_Time = col_annot_df$Treatment_Time,
  # Pick colors for each Treatment
  col = list(
    # Cancer = c("NSCLC" = "#114B5F",
    #                       "Breast Carcinoma" = "#1A936F",
    #                       "Esophagus Carcinoma" = "#88D498"
    # ),
    
    # Tissue = c("AU565 cancer cell line" = "#071E22", 
    #            "CTG-2534 PDX model" = "#1D7874", 
    #            "NCI-H1781 cancer cell line" = "#679289", 
    #            "NCI-H2170-YVMA CRISPR/Cas9 engineered cell line" = "#F4C095", 
    #            "OE19 cancer cell line" = "#EE2E31"
    # ),
    
    Drug = c(#"DMSO" = "#C9CBA3",
             "Natrosol" = "#FFE1A8",
             #"BI-1622" = "#E26D5C",
             "BI-4142" = "#723D46"#,
             #"Poziotinib" = "#472D30"
    ),
    Treatment_Time = c(
      #"0h" = "#FAF3DD",
      #"0.5h" = "#C8D5B9",
      #"2h" = "#8FC0A9",
      "6h" = "#68B0AB",
      #"16h" = "#4A7C59",
      "24h" = "#364958"
    )
  )
)

# ROW ANNOTATION
row_annot_file <- file.path(data_dir, "24_11.27_MT-SpaceBio_Pathways_Hm_Recent_meta_Arif_updated.csv")
row_annot_df <- readr::read_csv(row_annot_file)
row_annot_df <- row_annot_df %>%
  dplyr::filter(Pathway %in% rownames(gsva_results)) %>%
  tibble::column_to_rownames(var = "Pathway")

row_annot <- ComplexHeatmap::rowAnnotation(
  # Supply Comparison labels
  Pathway_Hierarchy = row_annot_df$Pathway_Hierarchy,
  # Pick colors for each column
  col = list(
    Pathway_Hierarchy = c(
      "Senescence"= "#a44e3f",                # rusty red (inflammaging / stress)
      "Innate Immune" = "#c24604",             # same as Immune anchor
      "Adaptive Immune" = "#d95d18",           # lighter burnt orange
      "mtDNA/dsRNA Innate Immune" = "#e1792a", # orange-brown
      "ISR" = "#705555",                       # stress gray/brown
      "UPR MT/ER" = "#b2786a",                 # pinkish brown (ER stress)
      "PANoptosis" = "#95342e",                # dark red (cell death)
      "OXPHOS" = "#5c755b",                    # olive green anchor
      "Metabolism" = "#9ba33b",                # green-yellow anchor
      "Antioxidant Defenses" = "#758d3d",      # earthy green
      "Cytosolic Protein Import" = "#6e7b4e",  # dull olive
      "Global Translation" = "#b5be4b",        # yellow-green
      "Cyto-Ribosome" = "#a6ac3e",             # similar to Metabolism
      "HIF" = "#4c5b52",                       # hypoxia dark green-gray
      "mTOR" = "#697b3c",                      # regulatory green
      "Mitophagy" = "#44553f",                 # dark moss green
      "RAAS" = "#914d1c"                       # brown-orange (vascular)
    )
  )
)

## Separate heatmaps for each cell line

# Generate col_split using the Comparison column
column_split <- factor(col_annot_df$Drug, levels = unique(col_annot_df$Drug))
row_split <- factor(row_annot_df$Pathway_Hierarchy, levels = unique(row_annot_df$Pathway_Hierarchy))

###########################

# Check if column annotation df and expression data are in the same order
all.equal(colnames(ctg_gsva_results), rownames(col_annot_df))
# Re-order data
ctg_gsva_results <- gsva_results %>%
  as.data.frame() %>%
  dplyr::select(rownames(col_annot_df)) %>%
  as.matrix()

# Create a color function based on standardized scale
color_func <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("navy", "white", "firebrick3")
)


library(ComplexHeatmap)
ctg_gsva_hm <- ComplexHeatmap::Heatmap(ctg_gsva_results, 
                                      name = "GSVA score",
                                      show_row_names = TRUE,
                                      row_names_side = "right",
                                      show_column_names = FALSE,
                                      cluster_rows = FALSE, 
                                      cluster_columns = FALSE,
                                      show_row_dend = FALSE,
                                      column_split = column_split,
                                      row_split = row_split,
                                      top_annotation = col_annot,
                                      left_annotation = row_annot,
                                      row_names_gp = gpar(fontsize = 10),  # Reduce font size for row names
                                      col = color_func#,
                                      #width = unit(20, "cm"),
                                      #height = unit(20, "cm")
)

ctg_gsva_heatmap <- ggplotify::as.ggplot(grid.grabExpr(draw(ctg_gsva_hm, 
                                                           heatmap_legend_side = "bottom", 
                                                           annotation_legend_side = "bottom", 
                                                           merge_legend = FALSE, 
                                                           padding = unit(c(5, 5, 5, 25), "mm")) # Add extra space for row names
))

# Save the Heatmap plot
ggsave(
  plot = ctg_gsva_heatmap,
  file.path(plots_dir, "NSCLC_CTG_CustomMitoSpacePathways_gsva_heatmap.png"),
  width = 8,
  height = 16
)

################################################################################
## NSLC NCI-H2170-YVMA CRISPR/Cas9 engineered cell line
# Set up column annotation from metadata
col_annot_df <- metadata %>%
  dplyr::filter(Cancer == "NSCLC") %>%
  dplyr::filter(Tissue == "NCI-H2170-YVMA CRISPR/Cas9 engineered cell line") %>%
  # Only select the columns of interest
  dplyr::select(Sample_ID, Tissue, Drug, Treatment_Time) %>%
  # Arrange by cluster and short_histology
  dplyr::arrange(Drug, Treatment_Time) %>%
  dplyr::arrange(factor(Treatment_Time, levels = c("0h", "0.5h", "2h", "16h"))) %>%
  dplyr::arrange(factor(Drug, levels = c("DMSO", "BI-1622", "Poziotinib"))) %>%
  # Store sample
  tibble::column_to_rownames("Sample_ID")


# Create the ComplexHeatmap column annotation object
col_annot <- ComplexHeatmap::HeatmapAnnotation(
  # Supply Treatment labels
  #Cancer = col_annot_df$Cancer,
  #Tissue = col_annot_df$Tissue,
  Drug = col_annot_df$Drug,
  Treatment_Time = col_annot_df$Treatment_Time,
  # Pick colors for each Treatment
  col = list(
    # Cancer = c("NSCLC" = "#114B5F",
    #                       "Breast Carcinoma" = "#1A936F",
    #                       "Esophagus Carcinoma" = "#88D498"
    # ),
    
    # Tissue = c("AU565 cancer cell line" = "#071E22", 
    #            "CTG-2534 PDX model" = "#1D7874", 
    #            "NCI-H1781 cancer cell line" = "#679289", 
    #            "NCI-H2170-YVMA CRISPR/Cas9 engineered cell line" = "#F4C095", 
    #            "OE19 cancer cell line" = "#EE2E31"
    # ),
    
    Drug = c("DMSO" = "#C9CBA3",
      #"Natrosol" = "#FFE1A8",
      "BI-1622" = "#E26D5C",
      #"BI-4142" = "#723D46",
      "Poziotinib" = "#472D30"
    ),
    Treatment_Time = c(
      "0h" = "#FAF3DD",
      "0.5h" = "#C8D5B9",
      "2h" = "#8FC0A9",
      #"6h" = "#68B0AB",
      "16h" = "#4A7C59"#,
      #"24h" = "#364958"
    )
  )
)

# ROW ANNOTATION
row_annot_file <- file.path(data_dir, "24_11.27_MT-SpaceBio_Pathways_Hm_Recent_meta_Arif_updated.csv")
row_annot_df <- readr::read_csv(row_annot_file)
row_annot_df <- row_annot_df %>%
  dplyr::filter(Pathway %in% rownames(gsva_results)) %>%
  tibble::column_to_rownames(var = "Pathway")

row_annot <- ComplexHeatmap::rowAnnotation(
  # Supply Comparison labels
  Pathway_Hierarchy = row_annot_df$Pathway_Hierarchy,
  # Pick colors for each column
  col = list(
    Pathway_Hierarchy = c(
      "Senescence"= "#a44e3f",                # rusty red (inflammaging / stress)
      "Innate Immune" = "#c24604",             # same as Immune anchor
      "Adaptive Immune" = "#d95d18",           # lighter burnt orange
      "mtDNA/dsRNA Innate Immune" = "#e1792a", # orange-brown
      "ISR" = "#705555",                       # stress gray/brown
      "UPR MT/ER" = "#b2786a",                 # pinkish brown (ER stress)
      "PANoptosis" = "#95342e",                # dark red (cell death)
      "OXPHOS" = "#5c755b",                    # olive green anchor
      "Metabolism" = "#9ba33b",                # green-yellow anchor
      "Antioxidant Defenses" = "#758d3d",      # earthy green
      "Cytosolic Protein Import" = "#6e7b4e",  # dull olive
      "Global Translation" = "#b5be4b",        # yellow-green
      "Cyto-Ribosome" = "#a6ac3e",             # similar to Metabolism
      "HIF" = "#4c5b52",                       # hypoxia dark green-gray
      "mTOR" = "#697b3c",                      # regulatory green
      "Mitophagy" = "#44553f",                 # dark moss green
      "RAAS" = "#914d1c"                       # brown-orange (vascular)
    )
  )
)

## Separate heatmaps for each cell line

# Generate col_split using the Comparison column
column_split <- factor(col_annot_df$Drug, levels = unique(col_annot_df$Drug))
row_split <- factor(row_annot_df$Pathway_Hierarchy, levels = unique(row_annot_df$Pathway_Hierarchy))

###########################

# Check if column annotation df and expression data are in the same order
all.equal(colnames(h2170_gsva_results), rownames(col_annot_df))
# Re-order data
h2170_gsva_results <- gsva_results %>%
  as.data.frame() %>%
  dplyr::select(rownames(col_annot_df)) %>%
  as.matrix()

# Create a color function based on standardized scale
color_func <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("navy", "white", "firebrick3")
)


library(ComplexHeatmap)
h2170_gsva_hm <- ComplexHeatmap::Heatmap(h2170_gsva_results, 
                                       name = "GSVA score",
                                       show_row_names = TRUE,
                                       row_names_side = "right",
                                       show_column_names = FALSE,
                                       cluster_rows = FALSE, 
                                       cluster_columns = FALSE,
                                       show_row_dend = FALSE,
                                       column_split = column_split,
                                       row_split = row_split,
                                       top_annotation = col_annot,
                                       left_annotation = row_annot,
                                       row_names_gp = gpar(fontsize = 10),  # Reduce font size for row names
                                       col = color_func#,
                                       #width = unit(20, "cm"),
                                       #height = unit(20, "cm")
)

h2170_gsva_heatmap <- ggplotify::as.ggplot(grid.grabExpr(draw(h2170_gsva_hm, 
                                                            heatmap_legend_side = "bottom", 
                                                            annotation_legend_side = "bottom", 
                                                            merge_legend = FALSE, 
                                                            padding = unit(c(5, 5, 5, 25), "mm")) # Add extra space for row names
))

# Save the Heatmap plot
ggsave(
  plot = h2170_gsva_heatmap,
  file.path(plots_dir, "NSCLC_H2170_CustomMitoSpacePathways_gsva_heatmap.png"),
  width = 8,
  height = 16
)

################################################################################
## NSLC NCI-H1781 cancer cell line
# Set up column annotation from metadata
col_annot_df <- metadata %>%
  dplyr::filter(Cancer == "NSCLC") %>%
  dplyr::filter(Tissue == "NCI-H1781 cancer cell line") %>%
  # Only select the columns of interest
  dplyr::select(Sample_ID, Tissue, Drug, Treatment_Time) %>%
  # Arrange by cluster and short_histology
  dplyr::arrange(Drug, Treatment_Time) %>%
  dplyr::arrange(factor(Treatment_Time, levels = c("6h", "24h"))) %>%
  dplyr::arrange(factor(Drug, levels = c("DMSO", "BI-4142", "Poziotinib"))) %>%
  # Store sample
  tibble::column_to_rownames("Sample_ID")


# Create the ComplexHeatmap column annotation object
col_annot <- ComplexHeatmap::HeatmapAnnotation(
  # Supply Treatment labels
  #Cancer = col_annot_df$Cancer,
  #Tissue = col_annot_df$Tissue,
  Drug = col_annot_df$Drug,
  Treatment_Time = col_annot_df$Treatment_Time,
  # Pick colors for each Treatment
  col = list(
    # Cancer = c("NSCLC" = "#114B5F",
    #                       "Breast Carcinoma" = "#1A936F",
    #                       "Esophagus Carcinoma" = "#88D498"
    # ),
    
    # Tissue = c("AU565 cancer cell line" = "#071E22", 
    #            "CTG-2534 PDX model" = "#1D7874", 
    #            "NCI-H1781 cancer cell line" = "#679289", 
    #            "NCI-H2170-YVMA CRISPR/Cas9 engineered cell line" = "#F4C095", 
    #            "OE19 cancer cell line" = "#EE2E31"
    # ),
    
    Drug = c("DMSO" = "#C9CBA3",
             #"Natrosol" = "#FFE1A8",
             #"BI-1622" = "#E26D5C",
             "BI-4142" = "#723D46",
             "Poziotinib" = "#472D30"
    ),
    Treatment_Time = c(
      #"0h" = "#FAF3DD",
      #"0.5h" = "#C8D5B9",
      #"2h" = "#8FC0A9",
      "6h" = "#68B0AB",
      #"16h" = "#4A7C59"#,
      "24h" = "#364958"
    )
  )
)

# ROW ANNOTATION
row_annot_file <- file.path(data_dir, "24_11.27_MT-SpaceBio_Pathways_Hm_Recent_meta_Arif_updated.csv")
row_annot_df <- readr::read_csv(row_annot_file)
row_annot_df <- row_annot_df %>%
  dplyr::filter(Pathway %in% rownames(gsva_results)) %>%
  tibble::column_to_rownames(var = "Pathway")

row_annot <- ComplexHeatmap::rowAnnotation(
  # Supply Comparison labels
  Pathway_Hierarchy = row_annot_df$Pathway_Hierarchy,
  # Pick colors for each column
  col = list(
    Pathway_Hierarchy = c(
      "Senescence"= "#a44e3f",                # rusty red (inflammaging / stress)
      "Innate Immune" = "#c24604",             # same as Immune anchor
      "Adaptive Immune" = "#d95d18",           # lighter burnt orange
      "mtDNA/dsRNA Innate Immune" = "#e1792a", # orange-brown
      "ISR" = "#705555",                       # stress gray/brown
      "UPR MT/ER" = "#b2786a",                 # pinkish brown (ER stress)
      "PANoptosis" = "#95342e",                # dark red (cell death)
      "OXPHOS" = "#5c755b",                    # olive green anchor
      "Metabolism" = "#9ba33b",                # green-yellow anchor
      "Antioxidant Defenses" = "#758d3d",      # earthy green
      "Cytosolic Protein Import" = "#6e7b4e",  # dull olive
      "Global Translation" = "#b5be4b",        # yellow-green
      "Cyto-Ribosome" = "#a6ac3e",             # similar to Metabolism
      "HIF" = "#4c5b52",                       # hypoxia dark green-gray
      "mTOR" = "#697b3c",                      # regulatory green
      "Mitophagy" = "#44553f",                 # dark moss green
      "RAAS" = "#914d1c"                       # brown-orange (vascular)
    )
  )
)

## Separate heatmaps for each cell line

# Generate col_split using the Comparison column
column_split <- factor(col_annot_df$Drug, levels = unique(col_annot_df$Drug))
row_split <- factor(row_annot_df$Pathway_Hierarchy, levels = unique(row_annot_df$Pathway_Hierarchy))

###########################

# Check if column annotation df and expression data are in the same order
all.equal(colnames(h1781_gsva_results), rownames(col_annot_df))
# Re-order data
h1781_gsva_results <- gsva_results %>%
  as.data.frame() %>%
  dplyr::select(rownames(col_annot_df)) %>%
  as.matrix()

# Create a color function based on standardized scale
color_func <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("navy", "white", "firebrick3")
)


library(ComplexHeatmap)
h1781_gsva_hm <- ComplexHeatmap::Heatmap(h1781_gsva_results, 
                                         name = "GSVA score",
                                         show_row_names = TRUE,
                                         row_names_side = "right",
                                         show_column_names = FALSE,
                                         cluster_rows = FALSE, 
                                         cluster_columns = FALSE,
                                         show_row_dend = FALSE,
                                         column_split = column_split,
                                         row_split = row_split,
                                         top_annotation = col_annot,
                                         left_annotation = row_annot,
                                         row_names_gp = gpar(fontsize = 10),  # Reduce font size for row names
                                         col = color_func#,
                                         #width = unit(20, "cm"),
                                         #height = unit(20, "cm")
)

h1781_gsva_heatmap <- ggplotify::as.ggplot(grid.grabExpr(draw(h1781_gsva_hm, 
                                                              heatmap_legend_side = "bottom", 
                                                              annotation_legend_side = "bottom", 
                                                              merge_legend = FALSE, 
                                                              padding = unit(c(5, 5, 5, 25), "mm")) # Add extra space for row names
))

# Save the Heatmap plot
ggsave(
  plot = h1781_gsva_heatmap,
  file.path(plots_dir, "NSCLC_H1781_CustomMitoSpacePathways_gsva_heatmap.png"),
  width = 8,
  height = 16
)

