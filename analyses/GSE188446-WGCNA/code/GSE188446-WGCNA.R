# Script to perform weighted correlation gene co-expression network analysis (WGCNA)
# Author: Shehbeel Arif
# Code adapted from ALSF (https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#4_Identifying_co-expression_gene_modules_with_WGCNA_-_RNA-seq)
# NASA GeneLab

# Load libraries
library(tidyverse)
library(DESeq2)
library(WGCNA)
library(ggpubr)
allowWGCNAThreads() # allow multi-threading (optional)

################################################################################
# SETTING UP DIRECTORIES & LOADING DATA
# Set up analysis directories
# Create the data folder if it doesn't exist
# if (!dir.exists("data")) {
#   dir.create("data")
# }

# # Define the file path to the plots directory
# plots_dir <- "plots"
# 
# # Create the plots folder if it doesn't exist
# if (!dir.exists(plots_dir)) {
#   dir.create(plots_dir)
# }
# 
# # Define the file path to the results directory
# results_dir <- "results"
# 
# # Create the results folder if it doesn't exist
# if (!dir.exists(results_dir)) {
#   dir.create(results_dir)
# }




## Set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "GSE188446-WGCNA")
data_dir <- file.path(root_dir, "data")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

# Specify dataset files
metadata_file <- file.path(data_dir, "GSE188446_meta.txt")
counts_file <- file.path(data_dir, "GSE188446_counts.txt")

# Load datasets
meta <- read.table(metadata_file, header=TRUE)
counts <- read.table(counts_file, header=TRUE, row.names=1)

################################################################################
## DATA PREPROCESSING

## Fix sample names in meta data
# Substitute "-" for "."
meta$Sample <- gsub("-", ".", meta$Sample)

# Rearrange counts column names by cluster data
counts <- counts[,meta$Sample]

# Convert the assigned cluster values into factors
meta$Treatment <- factor(meta$Treatment)

################################################################################
# QUALITY CONTROL - OUTLIER DETECTION
# Detect outlier genes
gsg <- goodSamplesGenes(t(counts))
summary(gsg)
gsg$allOK

table(gsg$goodGenes) # 24822 outlier genes detected
table(gsg$goodSamples) # No outlier samples detected

# Remove genes that are detectd as outliers
counts <- counts[gsg$goodGenes == TRUE,]

# Detect outlier samples using hierarchical clustering - method 1
htree <- hclust(dist(t(counts)), method = "average")
plot(htree)
# Sample outliers: 7316-496 (Ependymoma), 7316-230 (Ganglioglioma), 7316-88 (Ependymoma)

# Detect outlier samples using PCA - method 2

pca <- prcomp(t(counts))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# Exclude outlier samples
#samples.to.be.excluded <- c("7316-496")
#counts.subset <- counts[, !colnames(counts) %in% samples.to.be.excluded]

# Remove the outlier sample from colData as well
#clusters <- clusters[clusters$Sample_ID != "7316-496",]

################################################################################
# NORMALIZATION

# Normalize counts with DESeq
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~1) # not specifying the model

# Remove all genes with counts <15 in more than 75% of samples (0.75 * 14 samples = 10.5)
## Suggested by WGCNA on RNAseq FAQ
dds75 <- dds[rowSums(counts(dds) >= 15) > 11]
nrow(dds75) # 13498 genes left

# Perform variance stabilization
#dds_norm <- vst(dds75) # Did not work b/c it is less than 'nsub' rows
dds_norm <- varianceStabilizingTransformation(dds75)

# Transpose the counts
# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data

### Plot the normalized expressions
norm_counts_df <- data.frame(t(normalized_counts)) %>%
  mutate(
    Gene_id = row.names(t(normalized_counts))
  ) %>%
  pivot_longer(-Gene_id)

norm_counts_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "VST Normalized mRNA Expression",
    x = "Sample ID",
    y = "VST Normalized Expression"
  )

################################################################################
### WGCNA ANALYSIS ###
################################################################################
# Determine soft-threshold for scale-free network

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(normalized_counts,
                         #blockSize = 30,
                         powerVector = powers,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed", # genes are similar if they are strongly correlated; "signed" applied only to genes that are positively correlated
                         verbose = 5)

## Make plots of the soft-thresholds
par(mfrow = c(1,2));
cex1 = 0.8;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.80, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")


################################################################################
# Run WGCNA
chosen_power <- 14
temp_cor <- cor       
cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)

bwnet <- blockwiseModules(normalized_counts,
                          maxBlockSize = 5000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = chosen_power, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234 # there's some randomness associated with this calculation
                          # so we should set a seed
)


# Return cor function to original namespace
cor <- temp_cor

module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

# Look at number of genes per Modules
table(bwnet$colors)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21 
# 1059 3958 1043  801  494  312  310  292  284  265  221  213  212  186  182  176  175  162  131  111  110  110 
# 22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42   43 
# 104  103  103  101  100   97   89   89   85   81   81   74   72   70   70   68   67   66   62   62   62   62 
# 44   45   46   47   48   49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64   65 
# 56   56   54   51   50   50   47   47   41   40   39   39   39   36   35   32   32   31   29   27   25   24 
# 66   67 
# 22   21 

################################################################################
# WGCNA MODULE DENDROGRAM

# Convert labels to colors for plotting
mergedColors <- labels2colors(bwnet$colors)
unmergedColors <- labels2colors(bwnet$unmergedColors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet$dendrograms[[1]],
                    cbind(mergedColors[bwnet$blockGenes[[1]]], unmergedColors[bwnet$blockGenes[[1]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05 )


################################################################################
## WHICH MODULES HAVE BIGGEST DIFFERENCES ACROSS CLUSTER GROUPS?

# Check if SAMPLE_ID orders are the same as module_eigengenes
all.equal(meta$Sample, rownames(module_eigengenes))

# Create the design matrix from the `time_point` variable
des_mat <- model.matrix(~ meta$Treatment)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

# Save the Limma stats data to a CSV file
readr::write_csv(stats_df,
                 file = file.path(results_dir, "wgcna-limma.csv")
)

################################################################################
# CORRELATION BETWEEN MODULES AND CLUSTERS

# As a sanity check, let’s use ggplot to see what module 3’s eigengene looks like between treatment groups.

meta$sample_id <- gsub("-", ".", meta$sample_id)

module_df <- module_eigengenes %>%
  tibble::rownames_to_column("Sample") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(meta %>%
                      dplyr::select(sample_id, Treatment),
                    by = c("Sample" = "sample_id")
  )


# Reorder the Treatment Groups in module_df
module_df <- module_df %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Vehicle", "TAK", "TDM1", "Combo", "Resistance"))
  )

# Save the module eigengenes Sample_ID dataframe
readr::write_csv(module_df,
                 file = file.path(results_dir, "wgcna-sample-module.csv")
)

# Reorder modules so similar modules are next to each other
module_order <- names(module_df)

# Add clusters
temp_module_df <- module_df %>% dplyr::select("ME8":"Treatment")

# tidy & plot data
temp_module_df <- temp_module_df %>%
  pivot_longer(-Treatment) %>%
  mutate(
    module = factor(name, levels = module_order)
  )

temp_module_df %>% ggplot(., aes(x=Treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "#f7f7f7",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-Treatment Relationships", y = "Modules", x="Treatment", fill="corr")


################################################################################
# Make boxplot of Module Eigengene expression by clusters

# Create list of pairwise vectors to statistically compare clusters
my_comparisons <- list(c("Vehicle","TDM1"),
                       c("Vehicle","TAK"),
                       c("Vehicle","Combo"),
                       c("Vehicle","Resistance")
)

# Make plot
ME_plot <- 4
# Save plot
ggsave(file.path(plots_dir, "ME22_expression_in_treatments_boxplot.tiff"), ME_plot, width = 10, height = 6)


# ggplot(
#   module_df,
#   aes(
#     x = Treatment,
#     y = ME0,
#     color = Treatment
#   )
# ) +
#   labs(title = "Module Expression in Treatment Groups") +
#   # a boxplot with outlier points hidden (they will be in the sina plot)
#   geom_boxplot(#width = 0.2, 
#                outlier.shape = NA
#                ) +
#   # A sina plot to show all of the individual data points
#   ggforce::geom_sina(#maxwidth = 0.3
#                      ) +
#   # Add global p-value
#   stat_compare_means(#label.y = 0.35
#                      ) + # ME0=0.15, ME1=0.25, ME2=0.15, ME3=0.25
#   # Add pairwise comparisons p-values
#   #stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test') +
#   theme_classic()

################################################################################
## What genes are part of Module 62?
# Add Modules to genes
gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))
# Filter for genes part of only ME62
gene_module_key %>%
  dplyr::filter(module == "ME62")

# Save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path(results_dir, "wgcna_gene_to_module.tsv")
)


################################################################################
# MODULE-SPECIFIC HEATMAPS

# Make custom heatmap function
make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = meta,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its Sample_ID
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("Sample")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(Sample, Treatment) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "Sample") %>%
    # Arrange by patient and time point
    dplyr::arrange(Treatment) %>%
    # Store sample
    tibble::column_to_rownames("Sample")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    Treatment = col_annot_df$Treatment,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each Treatment group
    col = list(Treatment = c("Vehicle" = "#E64B35FF", "TAK" = "#4DBBD5FF", "TDM1" = "#00A087FF", "Combo" = "#3C5488FF", "Resistance" = "#F39B7FFF"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = TRUE,
                                     show_column_names = FALSE#,
                                     #column_title = "ME3 Gene Expression"
  )
  
  # Return heatmap
  return(heatmap)
}

# Make Module 3 heatmap
mod_3_heatmap <- make_module_heatmap(module_name = "ME14")

# Print out the plot
mod_3_heatmap


# Save this plot to PNG
png(file.path("results", "module_3_heatmap.png"))
mod_3_heatmap
#dev.off()

## For comparison look at other modules
#mod_2_heatmap <- make_module_heatmap(module_name = "ME2")
# Print out the plot
#mod_2_heatmap


################################################################################
## FINDING DRIVER GENES
chooseTopHubInEachModule(normalized_counts,
                         gene_module_key$module, 
                         power = chosen_power, 
                         type = "signed")

# ME0             ME1            ME10            ME11            ME12            ME13            ME14 
# "Ccnd1"         "Bud31"         "Cped1"         "Cops2"          "Gsap"         "Dram2"          "Mmp2" 
# ME15            ME16            ME17            ME18            ME19             ME2            ME20 
# "Sfmbt1"      "Ppp1r13l"         "Rad51"         "Wnt7b"          "Maea"        "Zfp652"      "Trappc6b" 
# ME21            ME22            ME23            ME24            ME25            ME26            ME27 
# "Irgm1"         "Wdr12"        "Klhl22"        "Dtnbp1"         "Ptpro"         "Igsf6"         "Smagp" 
# ME28            ME29             ME3            ME30            ME31            ME32            ME33 
# "Ldlrap1"         "Shmt1"         "Rnf44"         "Bspry"          "Krt5"          "Phf1"          "Fhit" 
# ME34            ME35            ME36            ME37            ME38            ME39             ME4 
# "Axin2"       "Slc12a6"         "Acsl4"         "Abhd2"        "Plxnb2"       "Slc12a9"         "Thoc2" 
# ME40            ME41            ME42            ME43            ME44            ME45            ME46 
# "Garnl3"        "Zfp788"      "Tmem184c"        "Sdcbp2"       "Gm28529"         "Wdfy4"      "Pisd-ps2" 
# ME47            ME48            ME49             ME5            ME50            ME51            ME52 
# "Gypc"         "Cldn3"         "Fbxw5"        "Map2k2"         "Mroh1"         "Mdfic"         "Cd164" 
# ME53            ME54            ME55            ME56            ME57            ME58            ME59 
# "Ccdc186" "1810010H24Rik"         "Prrg1"         "Nxpe5"       "Ccdc74a"        "Hdgfl3"        "Eif4a3" 
# ME6            ME60            ME61            ME62            ME63            ME64            ME65 
# "Fer1l4"        "Dnajb9"          "Chrd"         "Tinf2"         "Vps39"           "Kiz"       "Zfyve27" 
# ME66            ME67             ME7             ME8             ME9 
# "Acvrl1"         "Car5b"          "Osr1"       "Gm28729"         "Fbxw8" 


################################################################################
## EXAMINE EXPRESSION PROFILES

# pick out a few modules of interest here
modules_of_interest <- c("ME26", "ME27", "ME62")

# Pull out list of genes in that module
submod <- gene_module_key %>%
  subset(module %in% modules_of_interest)

row.names(gene_module_key) <- gene_module_key$gene

# Get normalized expression for those genes
temp_normalized_counts <- t(normalized_counts)
subexpr <- temp_normalized_counts[submod$gene,]

submod_df <- data.frame(subexpr) %>%
  mutate(
    gene = row.names(.)
  ) %>%
  pivot_longer(-gene) %>%
  mutate(
    module = gene_module_key[gene,]$module
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "Sample_ID",
       y = "Normalized Expression") #+
  #scale_color_manual(values=c("grey", "turquoise", "blue", "brown"))

################################################################################
## GENERATE NETWORKS
# The network file can be generated for Cytoscape or as an edge/vertices file.

genes_of_interest <- module_df %>%
  subset(colors %in% modules_of_interest)
submod

expr_of_interest <- normalized_counts[submod$gene_id,]


# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM <- TOMsimilarityFromExpr(t(subexpr),
                             power = chosen_power)

# Add gene names to row and columns
row.names(TOM) <- row.names(subexpr)
colnames(TOM) <- row.names(subexpr)

edge_list <- data.frame(TOM) %>%
  mutate(gene1 = row.names(.)) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  mutate(gene2 = gsub("*\\.","-",as.character(gene2))) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = submod[gene1,]$module,
    module2 = submod[gene2,]$module
  )


# Select rows if gene1 and gen2 are from same module
edge_list[edge_list$module1 == edge_list$module2,]
# Select rows if gene1 and gen2 have a correlation >0.1
new_edge_list <- edge_list[edge_list$correlation > 0.5,]

a <- new_edge_list$gene1
b <- new_edge_list$gene2
idx <- order(c(seq_along(a), seq_along(b)))
edges <- c(a,b)[idx]

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")

write_delim(as.data.frame(TOM),
            file = "TOM_edgelist.tsv",
            delim = "\t")

read.csv("edgelist.tsv", sep="\t")

# Creating a Network graph using igraph
# Tutorial: https://kateto.net/networks-r-igraph
library(igraph) # Load the igraph package
rm(list = ls()) # Remove all the objects we created so far.

g1 <- graph( edges=edges, directed=F) 

plot(g1)
