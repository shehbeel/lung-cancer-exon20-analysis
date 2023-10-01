# Differential expression analysis of GSE181673 dataset
# Author: Shehbeel Arif

# Load libraries
library(dplyr)
library(DESeq2)

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
## Ensembl Gene identifier conversion
#Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)

keytypes(org.Hs.eg.db)

# create a data frame that shows the mapped gene symbols along with the expression data for the respective Ensembl IDs
expression_df <- data.frame(
  gene_symbol = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = expression_df$GENE,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "SYMBOL",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
  # from the data frame
  dplyr::filter(!is.na(gene_symbol)) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(expression_df, by = c("Ensembl" = "GENE")) %>%
  dplyr::select(-c(Ensembl))

#########################################################################################
#######
## PREPROCESS THE DATA
expression_df <- expression_df %>%
  distinct(gene_symbol, .keep_all = TRUE) %>% # Remove duplicate genes
  tibble::column_to_rownames("gene_symbol")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$Sample_ID)

# Check if the expression and metadata are in the same order
all.equal(colnames(expression_df), metadata$Sample_ID)

# Make Treatment a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    Treatment = factor(Treatment, levels = c("CTG-2534_Natrosol_6h", 
                                             "CTG-2534_BI-4142_100mgkg_6h", 
                                             "CTG-2534_BI-4142_100mgkg_24h", 
                                             "NCI-H2170-YVMA_DMSO_0h", 
                                             "NCI-H2170-YVMA_Poziotinib_28nM_0.5h", 
                                             "NCI-H2170-YVMA_BI-1622_120nM_0.5h",
                                             "NCI-H2170-YVMA_Poziotinib_28nM_2h",
                                             "NCI-H2170-YVMA_BI-1622_120nM_2h",
                                             "NCI-H2170-YVMA_Poziotinib_28nM_16h",
                                             "NCI-H2170-YVMA_BI-1622_120nM_16h",
                                             "AU565_DMSO_6h",
                                             "AU565_Poziotinib_10nM_6h",
                                             "AU565_BI-4142_70nM_6h",
                                             "AU565_Poziotinib_10nM_24h",
                                             "AU565_BI-4142_70nM_24h",
                                             "NCI-H1781_DMSO_6h",
                                             "NCI-H1781_Poziotinib_30nM_6h",
                                             "NCI-H1781_BI-4142_300nM_6h",
                                             "NCI-H1781_Poziotinib_30nM_24h",
                                             "NCI-H1781_BI-4142_300nM_24h",
                                             "OE19_DMSO_6h",
                                             "OE19_Poziotinib_10nM_6h",
                                             "OE19_BI-4142_100nM_6h",
                                             "OE19_Poziotinib_10nM_24h",
                                             "OE19_BI-4142_100nM_24h"
                                             ))
  )

#########################################################################################

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

CTG_2534_BI_4142_6h_vs_Natrosol_results <- results(deseq_object, c("Treatment", "CTG-2534_BI-4142_100mgkg_6h", "CTG-2534_Natrosol_6h"))

# Sort and filter DESeq2 results table and convert to dataframe
CTG_2534_BI_4142_6h_vs_Natrosol_df <- CTG_2534_BI_4142_6h_vs_Natrosol_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))


# # Map ENSEMBL ID to GENE SYMBOL
# CTG_2534_BI_4142_6h_vs_Natrosol_df <- data.frame(
#   GENE = mapIds(
#     # Replace with annotation package for the organism relevant to your data
#     org.Hs.eg.db,
#     keys = CTG_2534_BI_4142_6h_vs_Natrosol_df$Gene,
#     # Replace with the type of gene identifiers in your data
#     keytype = "ENSEMBL",
#     # Replace with the type of gene identifiers you would like to map to
#     column = "SYMBOL",
#     # This will keep only the first mapped value for each Ensembl ID
#     multiVals = "first"
#   )
# ) %>%
#   # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
#   # from the data frame
#   dplyr::filter(!is.na(GENE)) %>%
#   # Make an `Ensembl` column to store the rownames
#   tibble::rownames_to_column("Ensembl") %>%
#   # Now let's join the rest of the expression data
#   dplyr::inner_join(CTG_2534_BI_4142_6h_vs_Natrosol_df, by = c("Ensembl" = "Gene")) %>%
#   dplyr::select(-c(Ensembl))

# Save results as CSV
readr::write_csv(
  CTG_2534_BI_4142_6h_vs_Natrosol_df,
  file.path(
    results_dir,
    "CTG_2534_BI_4142_6h_vs_Natrosol_DEG.csv"
  )
)

CTG_2534_BI_4142_24h_vs_Natrosol_results <- results(deseq_object, c("Treatment", "CTG-2534_BI-4142_100mgkg_24h", "CTG-2534_Natrosol_6h"))

# Sort and filter DESeq2 results table and convert to dataframe
CTG_2534_BI_4142_24h_vs_Natrosol_df <- CTG_2534_BI_4142_24h_vs_Natrosol_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  CTG_2534_BI_4142_24h_vs_Natrosol_df,
  file.path(
    results_dir,
    "CTG-2534_BI_4142_24h_vs_Natrosol_DEG.csv"
  )
)

CTG_2534_BI_4142_24h_vs_BI_4142_6h_results <- results(deseq_object, c("Treatment", "CTG-2534_BI-4142_100mgkg_24h", "CTG-2534_BI-4142_100mgkg_6h"))

# Sort and filter DESeq2 results table and convert to dataframe
CTG_2534_BI_4142_24h_vs_BI_4142_6h_df <- CTG_2534_BI_4142_24h_vs_BI_4142_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  CTG_2534_BI_4142_24h_vs_BI_4142_6h_df,
  file.path(
    results_dir,
    "CTG_2534_BI_4142_24h_vs_BI_4142_6h_DEG.csv"
  )
)

#########################################################################################

NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_results <- results(deseq_object, c("Treatment", "NCI-H2170-YVMA_Poziotinib_28nM_0.5h", "NCI-H2170-YVMA_DMSO_0h"))

NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_df <- NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_df,
  file.path(
    results_dir,
    "NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_DEG.csv"
  )
)

NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_results <- results(deseq_object, c("Treatment", "NCI-H2170-YVMA_Poziotinib_28nM_16h", "NCI-H2170-YVMA_DMSO_0h"))

NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_df <- NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_df,
  file.path(
    results_dir,
    "NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_DEG.csv"
  )
)

NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_results <- results(deseq_object, c("Treatment", "NCI-H2170-YVMA_Poziotinib_28nM_16h", "NCI-H2170-YVMA_Poziotinib_28nM_0.5h"))

NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_df <- NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_df,
  file.path(
    results_dir,
    "NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_DEG.csv"
  )
)

#########################################################################################

NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_results <- results(deseq_object, c("Treatment", "NCI-H2170-YVMA_BI-1622_120nM_0.5h", "NCI-H2170-YVMA_DMSO_0h"))

NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_df <- NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_df,
  file.path(
    results_dir,
    "NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_DEG.csv"
  )
)


NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_results <- results(deseq_object, c("Treatment", "NCI-H2170-YVMA_BI-1622_120nM_16h", "NCI-H2170-YVMA_DMSO_0h"))

NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_df <- NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_df,
  file.path(
    results_dir,
    "NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_DEG.csv"
  )
)

NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_results <- results(deseq_object, c("Treatment", "NCI-H2170-YVMA_BI-1622_120nM_16h", "NCI-H2170-YVMA_BI-1622_120nM_0.5h"))

NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_df <- NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_df,
  file.path(
    results_dir,
    "NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_DEG.csv"
  )
)

#########################################################################################

AU565_Poziotinib_10nM_6h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "AU565_Poziotinib_10nM_6h", "AU565_DMSO_6h"))

AU565_Poziotinib_10nM_6h_vs_DMSO_6h_df <- AU565_Poziotinib_10nM_6h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  AU565_Poziotinib_10nM_6h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "AU565_Poziotinib_10nM_6h_vs_DMSO_6h_DEG.csv"
  )
)

AU565_Poziotinib_10nM_24h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "AU565_Poziotinib_10nM_24h", "AU565_DMSO_6h"))

AU565_Poziotinib_10nM_24h_vs_DMSO_6h_df <- AU565_Poziotinib_10nM_24h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  AU565_Poziotinib_10nM_24h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "AU565_Poziotinib_10nM_24h_vs_DMSO_6h_DEG.csv"
  )
)

AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_results <- results(deseq_object, c("Treatment", "AU565_Poziotinib_10nM_24h", "AU565_Poziotinib_10nM_6h"))

AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_df <- AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_df,
  file.path(
    results_dir,
    "AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_DEG.csv"
  )
)

#########################################################################################

AU565_BI_4142_70nM_6h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "AU565_BI-4142_70nM_6h", "AU565_DMSO_6h"))

AU565_BI_4142_70nM_6h_vs_DMSO_6h_df <- AU565_BI_4142_70nM_6h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  AU565_BI_4142_70nM_6h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "AU565_BI_4142_70nM_6h_vs_DMSO_6h_DEG.csv"
  )
)


AU565_BI_4142_70nM_24h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "AU565_BI-4142_70nM_24h", "AU565_DMSO_6h"))

AU565_BI_4142_70nM_24h_vs_DMSO_6h_df <- AU565_BI_4142_70nM_24h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  AU565_BI_4142_70nM_24h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "AU565_BI_4142_70nM_24h_vs_DMSO_6h_DEG.csv"
  )
)

AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_results <- results(deseq_object, c("Treatment", "AU565_BI-4142_70nM_24h", "AU565_BI-4142_70nM_6h"))

AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_df <- AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_df,
  file.path(
    results_dir,
    "AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_DEG.csv"
  )
)
#########################################################################################

NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "NCI-H1781_Poziotinib_30nM_6h", "NCI-H1781_DMSO_6h"))

NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_df <- NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_DEG.csv"
  )
)

NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "NCI-H1781_Poziotinib_30nM_24h", "NCI-H1781_DMSO_6h"))

NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_df <- NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_DEG.csv"
  )
)

NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_results <- results(deseq_object, c("Treatment", "NCI-H1781_Poziotinib_30nM_24h", "NCI-H1781_Poziotinib_30nM_6h"))

NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_df <- NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_df,
  file.path(
    results_dir,
    "NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_DEG.csv"
  )
)

#########################################################################################
NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "NCI-H1781_BI-4142_300nM_6h", "NCI-H1781_DMSO_6h"))

NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_df <- NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_DEG.csv"
  )
)

NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "NCI-H1781_BI-4142_300nM_24h", "NCI-H1781_DMSO_6h"))

NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_df <- NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_DEG.csv"
  )
)

NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_results <- results(deseq_object, c("Treatment", "NCI-H1781_BI-4142_300nM_24h", "NCI-H1781_BI-4142_300nM_6h"))

NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_df <- NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_df,
  file.path(
    results_dir,
    "NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_DEG.csv"
  )
)

#########################################################################################
OE19_Poziotinib_10nM_6h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "OE19_Poziotinib_10nM_6h", "OE19_DMSO_6h"))

OE19_Poziotinib_10nM_6h_vs_DMSO_6h_df <- OE19_Poziotinib_10nM_6h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  OE19_Poziotinib_10nM_6h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "OE19_Poziotinib_10nM_6h_vs_DMSO_6h_DEG.csv"
  )
)

OE19_Poziotinib_10nM_24h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "OE19_Poziotinib_10nM_24h", "OE19_DMSO_6h"))

OE19_Poziotinib_10nM_24h_vs_DMSO_6h_df <- OE19_Poziotinib_10nM_24h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  OE19_Poziotinib_10nM_24h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "OE19_Poziotinib_10nM_24h_vs_DMSO_6h_DEG.csv"
  )
)

OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_results <- results(deseq_object, c("Treatment", "OE19_Poziotinib_10nM_24h", "OE19_Poziotinib_10nM_6h"))

OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_df <- OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_df,
  file.path(
    results_dir,
    "OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_DEG.csv"
  )
)

#########################################################################################
OE19_BI_4142_100nM_6h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "OE19_BI-4142_100nM_6h", "OE19_DMSO_6h"))

OE19_BI_4142_100nM_6h_vs_DMSO_6h_df <- OE19_BI_4142_100nM_6h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  OE19_BI_4142_100nM_6h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "OE19_BI_4142_100nM_6h_vs_DMSO_6h_DEG.csv"
  )
)

OE19_BI_4142_100nM_24h_vs_DMSO_6h_results <- results(deseq_object, c("Treatment", "OE19_BI-4142_100nM_24h", "OE19_DMSO_6h"))

OE19_BI_4142_100nM_24h_vs_DMSO_6h_df <- OE19_BI_4142_100nM_24h_vs_DMSO_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  OE19_BI_4142_100nM_24h_vs_DMSO_6h_df,
  file.path(
    results_dir,
    "OE19_BI_4142_100nM_24h_vs_DMSO_6h_DEG.csv"
  )
)

OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_results <- results(deseq_object, c("Treatment", "OE19_BI-4142_100nM_24h", "OE19_BI-4142_100nM_6h"))

OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_df <- OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save results as CSV
readr::write_csv(
  OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_df,
  file.path(
    results_dir,
    "OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_DEG.csv"
  )
)

################################################################################
## Create volcano plot
## CTG_2534_BI_4142_6h_vs_Natrosol
CTG_2534_BI_4142_6h_vs_Natrosol_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  CTG_2534_BI_4142_6h_vs_Natrosol_df,
  lab = CTG_2534_BI_4142_6h_vs_Natrosol_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
CTG_2534_BI_4142_6h_vs_Natrosol_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = CTG_2534_BI_4142_6h_vs_Natrosol_volcano_plot,
  file.path(plots_dir, "CTG_2534_BI_4142_6h_vs_Natrosol_volcano_plot.png"),
  scale=1.5
)

## CTG_2534_BI_4142_24h_vs_Natrosol
CTG_2534_BI_4142_24h_vs_Natrosol_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  CTG_2534_BI_4142_24h_vs_Natrosol_df,
  lab = CTG_2534_BI_4142_24h_vs_Natrosol_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
CTG_2534_BI_4142_24h_vs_Natrosol_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = CTG_2534_BI_4142_24h_vs_Natrosol_volcano_plot,
  file.path(plots_dir, "CTG_2534_BI_4142_24h_vs_Natrosol_volcano_plot.png"),
  scale=1.5
)

## CTG_2534_BI_4142_24h_vs_BI_4142_6h
CTG_2534_BI_4142_24h_vs_BI_4142_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  CTG_2534_BI_4142_24h_vs_BI_4142_6h_df,
  lab = CTG_2534_BI_4142_24h_vs_BI_4142_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
CTG_2534_BI_4142_24h_vs_BI_4142_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = CTG_2534_BI_4142_24h_vs_BI_4142_6h_volcano_plot,
  file.path(plots_dir, "CTG_2534_BI_4142_24h_vs_BI_4142_6h_volcano_plot.png"),
  scale=1.5
)

## NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h
NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_df,
  lab = NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_volcano_plot,
  file.path(plots_dir, "NCI_H2170_YVMA_Poziotinib_28nM_0.5h_vs_DMSO_0h_volcano_plot.png"),
  scale=1.5
)

## NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h
NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_df,
  lab = NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_volcano_plot,
  file.path(plots_dir, "NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_DMSO_0h_volcano_plot.png"),
  scale=1.5
)

## NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h
NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_df,
  lab = NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_volcano_plot,
  file.path(plots_dir, "NCI_H2170_YVMA_Poziotinib_28nM_16h_vs_Poziotinib_28nM_0.5h_volcano_plot.png"),
  scale=1.5
)

## NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h
NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_df,
  lab = NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_volcano_plot,
  file.path(plots_dir, "NCI_H2170_YVMA_BI_1622_120nM_0.5h_vs_DMSO_0h_volcano_plot.png"),
  scale=1.5
)

## NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h
NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_df,
  lab = NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_volcano_plot,
  file.path(plots_dir, "NCI_H2170_BI_1622_120nM_16h_vs_YVMA_DMSO_0h_volcano_plot.png"),
  scale=1.5
)

## NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h
NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_df,
  lab = NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_volcano_plot,
  file.path(plots_dir, "NCI_H2170_YVMA_BI_1622_120nM_16h_vs_BI_1622_120nM_0.5h_volcano_plot.png"),
  scale=1.5
)

## AU565_Poziotinib_10nM_6h_vs_DMSO_6h
AU565_Poziotinib_10nM_6h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  AU565_Poziotinib_10nM_6h_vs_DMSO_6h_df,
  lab = AU565_Poziotinib_10nM_6h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
AU565_Poziotinib_10nM_6h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = AU565_Poziotinib_10nM_6h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "AU565_Poziotinib_10nM_6h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## AU565_Poziotinib_10nM_24h_vs_DMSO_6h
AU565_Poziotinib_10nM_24h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  AU565_Poziotinib_10nM_24h_vs_DMSO_6h_df,
  lab = AU565_Poziotinib_10nM_24h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
AU565_Poziotinib_10nM_24h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = AU565_Poziotinib_10nM_24h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "AU565_Poziotinib_10nM_24h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h
AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_df,
  lab = AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_volcano_plot,
  file.path(plots_dir, "AU565_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_volcano_plot.png"),
  scale=1.5
)

## AU565_BI_4142_70nM_6h_vs_DMSO_6h
AU565_BI_4142_70nM_6h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  AU565_BI_4142_70nM_6h_vs_DMSO_6h_df,
  lab = AU565_BI_4142_70nM_6h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
AU565_BI_4142_70nM_6h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = AU565_BI_4142_70nM_6h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "AU565_BI_4142_70nM_6h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## AU565_BI_4142_70nM_24h_vs_DMSO_6h
AU565_BI_4142_70nM_24h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  AU565_BI_4142_70nM_24h_vs_DMSO_6h_df,
  lab = AU565_BI_4142_70nM_24h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
AU565_BI_4142_70nM_24h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = AU565_BI_4142_70nM_24h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "AU565_BI_4142_70nM_24h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h
AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_df,
  lab = AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_volcano_plot,
  file.path(plots_dir, "AU565_BI_4142_70nM_24h_vs_BI_4142_70nM_6h_volcano_plot.png"),
  scale=1.5
)

## NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h
NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_df,
  lab = NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "NCI_H1781_Poziotinib_30nM_6h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h
NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_df,
  lab = NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "NCI_H1781_Poziotinib_30nM_24h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h
NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_df,
  lab = NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_volcano_plot,
  file.path(plots_dir, "NCI_H1781_Poziotinib_30nM_24h_vs_Poziotinib_30nM_6h_volcano_plot.png"),
  scale=1.5
)

## NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h
NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_df,
  lab = NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "NCI_H1781_BI_4142_300nM_6h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h
NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_df,
  lab = NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "NCI_H1781_BI_4142_300nM_24h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h
NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_df,
  lab = NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_volcano_plot,
  file.path(plots_dir, "NCI_H1781_BI_4142_300nM_24h_vs_BI_4142_300nM_6h_volcano_plot.png"),
  scale=1.5
)

## OE19_Poziotinib_10nM_6h_vs_DMSO_6h
OE19_Poziotinib_10nM_6h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  OE19_Poziotinib_10nM_6h_vs_DMSO_6h_df,
  lab = OE19_Poziotinib_10nM_6h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
OE19_Poziotinib_10nM_6h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = OE19_Poziotinib_10nM_6h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "OE19_Poziotinib_10nM_6h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)


## OE19_Poziotinib_10nM_24h_vs_DMSO_6h
OE19_Poziotinib_10nM_24h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  OE19_Poziotinib_10nM_24h_vs_DMSO_6h_df,
  lab = OE19_Poziotinib_10nM_24h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
OE19_Poziotinib_10nM_24h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = OE19_Poziotinib_10nM_24h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "OE19_Poziotinib_10nM_24h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h
OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_df,
  lab = OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_volcano_plot,
  file.path(plots_dir, "OE19_Poziotinib_10nM_24h_vs_Poziotinib_10nM_6h_volcano_plot.png"),
  scale=1.5
)

## OE19_BI_4142_100nM_6h_vs_DMSO_6h
OE19_BI_4142_100nM_6h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  OE19_BI_4142_100nM_6h_vs_DMSO_6h_df,
  lab = OE19_BI_4142_100nM_6h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
OE19_BI_4142_100nM_6h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = OE19_BI_4142_100nM_6h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "OE19_BI_4142_100nM_6h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## OE19_BI_4142_100nM_24h_vs_DMSO_6h
OE19_BI_4142_100nM_24h_vs_DMSO_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  OE19_BI_4142_100nM_24h_vs_DMSO_6h_df,
  lab = OE19_BI_4142_100nM_24h_vs_DMSO_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
OE19_BI_4142_100nM_24h_vs_DMSO_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = OE19_BI_4142_100nM_24h_vs_DMSO_6h_volcano_plot,
  file.path(plots_dir, "OE19_BI_4142_100nM_24h_vs_DMSO_6h_volcano_plot.png"),
  scale=1.5
)

## OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h
OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_df,
  lab = OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_volcano_plot,
  file.path(plots_dir, "OE19_BI_4142_100nM_24h_vs_BI_4142_100nM_6h_volcano_plot.png"),
  scale=1.5
)


