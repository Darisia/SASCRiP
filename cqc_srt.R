

# required packages
library(tidyverse)
library(Seurat)
library(patchwork)
library(gridExtra)
library(Matrix)
library(stats)

#################################################################################
# General R function to load an RData file using a specified object name
#################################################################################
loadRData <- function(filename) {
  load(filename)
  get(ls()[ls() != "filename"])
}
####################################################################################
# R function to replace the ENSEMBL IDs with hgnc gene names within a Seurat object
####################################################################################
ensg_to_hgnc <- function(
  ensg_gname_path,
  read_in_seurat_object,
  sample_ID
){
  print("Seurat object contains ENSEMBL gene names")
  # Set the path to the ensg_gname.tsv file
  print("Looking for tsv file containing ENSG names and HGNC gene symbols")
  # read in the ensg_gname.tsv file
  ensg_gname <- read_tsv(ensg_gname_path, col_names = FALSE)
  # Add colnames to the ensg_gname file
  colnames(ensg_gname) <- c("ENSEMBL_ID", "gene_name")
  # Assign read_in_seurat_object to sample_ID_srt
  sample_ID_srt <- read_in_seurat_object
  # Update Seurat object if saved as an old seurat object
  sample_ID_srt <- UpdateSeuratObject(sample_ID_srt)
  # Generate the dataframe that will be used to match and replace the ENSEMBL IDs
  pmatch_output <- data.frame(pmatch(ensg_gname$ENSEMBL_ID, rownames(sample_ID_srt@assays[["RNA"]]@counts)))
  colnames(pmatch_output) <- "gene_match"
  pmatch_output$match_index <- rownames(pmatch_output)
  pmatch_output <- pmatch_output[order(pmatch_output$gene_match), ]
  pmatch_output <- filter(pmatch_output, !is.na(gene_match))
  # Match and replace the ENSEMBL IDs in the relevant slots in the seurat object with the hgnc gene symbol
  print("Replacing ENSG names with corresponding HGNC gene symbols")
  rownames(sample_ID_srt@assays[["RNA"]]@counts) <- ensg_gname$gene_name[as.numeric(pmatch_output$match_index)]
  rownames(sample_ID_srt@assays[["RNA"]]@data) <- ensg_gname$gene_name[as.numeric(pmatch_output$match_index)]
  # Return the Seurat object with hgnc symbols as rownames
  return(sample_ID_srt)
}
######################################################
# Read in and check the seurat object
######################################################
read_check_srt <- function(
  seurat_object,
  sample_ID,
  ensg_gname_path
){
  # Read in the saved seurat object
  if (endsWith(seurat_object, ".rds") == TRUE | endsWith(seurat_object, "RDS") == TRUE) {
    sample_ID_srt <- readRDS(file = seurat_object)
  } else {
    sample_ID_srt <- loadRData(seurat_object)
  }
  # Update Seurat object if saved as an old seurat object
  sample_ID_srt <- UpdateSeuratObject(sample_ID_srt)
  # Check if rownames use hgnc symbols or ENSEMBL IDs - change to hgnc symbol if ENSEMBL IDs are used
  if (startsWith(rownames(sample_ID_srt@assays[["RNA"]]@counts)[1], "ENSG") == TRUE) {
    sample_ID_srt <- ensg_to_hgnc(ensg_gname_path, sample_ID_srt, sample_ID)
  }
  # Remove cell barcodes with counts of 0
  sample_ID_srt <- subset(sample_ID_srt, subset = nCount_RNA > 0)
  # Check if the mito.percent column is in the meta.data slot of the Seurat object
  if ("mito.percent" %in% colnames(sample_ID_srt@meta.data)) {
    sample_ID_srt <- sample_ID_srt
  } else {
    sample_ID_srt[['mito.percent']] <- PercentageFeatureSet(sample_ID_srt, pattern = '^MT-')
  }
  return(sample_ID_srt)
}
######################################################################
# Convert information from the Seurat object into a mtx count matrix
######################################################################
srt_to_mtx <- function(
  seurat_object,
  sample_ID,
  output_folder
){
  print("Converting Seurat Object to an mtx matrix")
  sample_ID_srt <- seurat_object
  sample_ID_sparse <- as.sparse(GetAssayData(sample_ID_srt, slot = 'counts'))
  writeMM(sample_ID_sparse, sprintf('%s/%s_QC_matrix_output/matrix.mtx', output_folder, sample_ID))
  sample_ID_sparse_rows <- rownames(sample_ID_sparse)
  sample_ID_sparse_cols <- colnames(sample_ID_sparse)
  write.table(sample_ID_sparse_rows, file = sprintf('%s/%s_QC_matrix_output/genes.tsv', output_folder, sample_ID), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(sample_ID_sparse_cols, file = sprintf('%s/%s_QC_matrix_output/barcodes.tsv', output_folder, sample_ID), row.names = FALSE, col.names = FALSE, quote = FALSE)
}
#################################################
# Determine the cutoff for total genes detected
#################################################
out_detect <- function(
  seurat_object,
  sample_ID,
  gene_lower,
  mitochondria_percent,
  outlier_method,
  nMADs,
  nSD
) {
  sample_ID_srt <- seurat_object
  # Remove poor quality cells and save the remaining information in a dataframe
  sample_ID_mtdata <- sample_ID_srt@meta.data
  sample_ID_mtdata <- filter(sample_ID_mtdata, nFeature_RNA > as.numeric(gene_lower) | mito.percent < as.numeric(mitochondria_percent))
  # Calculate outlier cut off based on outlier_method
  if (outlier_method == "MAD") {
    mad_value <- mad(sample_ID_mtdata$nFeature_RNA)
    n_mads_value <- mad_value * as.numeric(nMADs)
    mads_from_median <- median(sample_ID_mtdata$nFeature_RNA) + n_mads_value
    higher_outlier_value <- mads_from_median
  } else if (outlier_method == "SD") {
    sd_value <- sd(sample_ID_mtdata$nFeature_RNA)
    n_sd_value <- sd_value * as.numeric(nSD)
    sd_from_mean <- mean(sample_ID_mtdata$nFeature_RNA) + n_sd_value
    higher_outlier_value <- sd_from_mean
  }
  return(higher_outlier_value)
}
####################################################################
# Generate a usable Seurat object containing all major QC metrics
####################################################################
gen_srt <- function(
  input_counts,
  sample_ID,
  output_folder,
  ensg_gname_path,
  gene_column
){
  print("Generating the Seurat object from a gene-count matrix")
  # import the data according to the file type
  if (endsWith(input_counts, ".h5") == TRUE) {
    sample_ID_data <- Read10X_h5(filename = input_counts)
  } else {
    sample_ID_data <- Read10X(data.dir = input_counts,
                              gene.column = as.double(gene_column))
  }
  # Convert the data into a Seurat Object
  sample_ID_srt <- CreateSeuratObject(counts = sample_ID_data,
                                  project = sample_ID)
  # Remove barcodes that are not attached to any cells
  sample_ID_srt <- subset(sample_ID_srt, subset = nCount_RNA > 0)

  # Check the gene names of the Seurat object
  if (startsWith(rownames(sample_ID_srt@assays[["RNA"]]@counts)[1], "ENSG") == TRUE) {
    sample_ID_srt <- ensg_to_hgnc(ensg_gname_path, sample_ID_srt, sample_ID)
  }
  # Calculate the mitochondrial percentage per cell
  print("Calculating fraction of mitochondrial genes using '^MT-'")
  sample_ID_srt[['mito.percent']] <- PercentageFeatureSet(sample_ID_srt, pattern = '^MT-')
  # Save the preQC seurat object to a file using saveRDS
  print("Saving Seurat object in .rds format to specified directory")
  output_preQC_seurat <- sprintf("%s/%s_preQC_seurat.rds", output_folder, sample_ID)
  saveRDS(
    sample_ID_srt,
    file = output_preQC_seurat
  )
  return(output_preQC_seurat)
}
###################################################################
# Generate the QC violin plots from data in a preQC seurat object
###################################################################
vln_srt <- function(
  seurat_object,
  sample_ID,
  output_folder,
  ensg_gname_path
){

  # Read in the saved seurat object
  sample_ID_srt <- read_check_srt(seurat_object, sample_ID, ensg_gname_path)

  # Produce a violin plot of the QC metrics
  sample_ID_vln_feature <- VlnPlot(sample_ID_srt, features = 'nFeature_RNA') +
    scale_fill_manual(values = 'steelblue2') +
    labs(title = '',
         y = 'Number of unique genes detected',
         x = '',
         subtitle = '') +
    theme(axis.ticks.length.x = unit(0, 'cm'),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 16,
                                      margin = margin(l = 1, unit = 'cm'), vjust = 5),
          axis.text.y = element_text(size = 14),
          legend.position = 'none')

  sample_ID_vln_UMI <- VlnPlot(sample_ID_srt, features = 'nCount_RNA') +
    scale_fill_manual(values = 'steelblue2') +
    labs(title = '',
         y = 'Number of UMIs detected',
         x = '',
         subtitle = '') +
    theme(axis.ticks.length.x = unit(0, 'cm'),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 16,
                                      margin = margin(l = 1, unit = 'cm'), vjust = 5),
          axis.text.y = element_text(size = 14),
          legend.position = 'none')

  sample_ID_vln_mt <- VlnPlot(sample_ID_srt, features = 'mito.percent') +
    scale_fill_manual(values = 'steelblue2') +
    labs(title = '',
         y = 'Fraction of mitochondrial genes (%)',
         x = '',
         subtitle = '') +
    theme(axis.ticks.length.x = unit(0, 'cm'),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 16,
                                      margin = margin(l = 1, unit = 'cm'), vjust = 5),
          axis.text.y = element_text(size = 14),
          legend.position = 'none')

  sample_ID_vln_metrics <- grid.arrange(sample_ID_vln_feature, sample_ID_vln_UMI, sample_ID_vln_mt,
                                        ncol = 3)
  # Save the violin plot to a specified output directory
  ggsave(sprintf('%s/%s_vln_metrics.png',
                 output_folder,
                 sample_ID),
         sample_ID_vln_metrics,
         width = 15,
         height = 8,
         dpi = 600)
}
############################################################
# Generate the combined QC plot containing all QC metrics
############################################################
combined_qc <- function(
  seurat_object,
  sample_ID,
  output_folder,
  ensg_gname_path
){

  # Read in and check the input seurat object
  sample_ID_srt <- read_check_srt(seurat_object, sample_ID, ensg_gname_path)
  # Use the meta.data dataframe from the Seurat object to make the combined QC plot
  combined_met_plot <- ggplot(sample_ID_srt@meta.data,
                              aes(x = nFeature_RNA,
                                  y = nCount_RNA,
                                  colour = mito.percent)) +
    geom_point() +
    scale_color_gradient(low = 'blue', high = 'red', guide = 'colourbar') +
    labs(x = 'Number of unique genes detected',
         y = 'Number of UMIs detected',
         colour = 'Mitochondrial genes (%)') +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 20),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 20),
          axis.text.x = element_text(size = 14,
                                     colour = 'black'),
          axis.text.y = element_text(size = 14,
                                     colour = 'black'),
          legend.text = element_text(size = 12,
                                     colour = 'black'),
          legend.title = element_text(size = 14,
                                      colour = 'black',
                                      angle = 90,
                                      margin = margin(t = -0.9, unit = 'cm')),
          legend.text.align = 0.5,
          axis.line = element_line(colour = 'black'),
          strip.background = element_rect(fill = 'white', colour = 'black'),
          strip.text = element_text(size = 14)
    ) +
    guides(colour = guide_colorbar(title.position = 'right'))

  # Save the combined QC plot to a specified output directory
  ggsave(sprintf('%s/%s_all_cell_metrics.png',
                 output_folder,
                 sample_ID),
         combined_met_plot,
         width = 15,
         height = 10,
         dpi = 600)
}

##############################################################################
# Lower gene count and mitochondiral percentage scatter plots
##############################################################################

damaged_scatter <- function(
  seurat_object,
  sample_ID,
  output_folder,
  ensg_gname_path,
  gene_lower,
  mitochondria_percent
) {

  # Read in and check the input seurat object
  sample_ID_srt <- read_check_srt(seurat_object, sample_ID, ensg_gname_path)

  # Extract plotting information and create data frame
  sample_ID_preQC_df <- data.frame(sample_ID_srt@meta.data$nFeature_RNA,
                               sample_ID_srt@meta.data$nCount_RNA,
                               sample_ID_srt@meta.data$mito.percent)

  # Rename the columns in the preQC dataframe
  colnames(sample_ID_preQC_df) <- c("gene_count", "UMI_count", "mito_percent")

  # Add another column deciding whether to keep or Remove cells based on mt percentage
  mitochondria_percent <- as.double(mitochondria_percent)
  sample_ID_preQC_df$mito_filt <- ifelse(sample_ID_preQC_df$mito_percent>mitochondria_percent, "Remove", "Keep")

  # Arrange the mito_filt column where Keep cells are first and Remove cells are after
  sample_ID_preQC_df <- arrange(sample_ID_preQC_df, mito_filt)

  # Create the plot to flag gene counts below 200
  gene_lower <- as.double(gene_lower)
  sample_ID_preQC_scat <- ggplot(sample_ID_preQC_df, aes(x = gene_count,
                                                   y = UMI_count)) +
  geom_point(aes(colour = ifelse(gene_count <= gene_lower, "Remove", "Keep" ))) +
  scale_colour_manual(values = c("blue", "red"),
                      labels = c("Healthy", "Low-quality"))

  sample_ID_preQC_scat <- sample_ID_preQC_scat +
    labs(x = "Total number of unique genes",
       y = "Total number of UMIs detected",
       colour = "") +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                    vjust = -5,
                                    size = 16),
               axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                           vjust = 5,
                                           size = 16),
               axis.text.x = element_text(size = 12,
                                          colour = 'black'),
               axis.text.y = element_text(size = 12,
                                          colour = 'black'),
               legend.text = element_text(size = 12,
                                          colour = 'black'),
               legend.title = element_text(size = 14,
                                           colour = 'black'),
               legend.text.align = 0.5,
               plot.margin = unit(c(1,1,1,1), "cm")) +
    theme(panel.background = element_blank(),
         axis.line = element_line(colour = "black"))

  ggsave(sprintf("%s/%s_lower_gene.png",
                output_folder,
                sample_ID),
        sample_ID_preQC_scat,
        width = 210,
        height = 150,
        units = "mm",
        dpi = 900)

  # Create the plot for the mitochondrial genes percentage after lower genes
  sample_ID_preQC_mito <- sample_ID_preQC_df %>%
  filter(gene_count > gene_lower) %>%
    ggplot(aes(x = gene_count,
               y = mito_percent,
               colour = mito_filt)) +
   geom_point() +
   scale_colour_manual(values = c("blue", "#ff7f00"),
                        labels = c("Healthy", "Low-quality")) +
  labs(x = "Total number of unique genes detected",
       y = "Percentage of mitochondrial genes detected",
       colour = "") +
  theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                       vjust = -5,
                                       size = 16),
                  axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                              vjust = 5,
                                              size = 16),
                  axis.text.x = element_text(size = 12,
                                             colour = 'black'),
                  axis.text.y = element_text(size = 12,
                                             colour = 'black'),
                  legend.text = element_text(size = 12,
                                             colour = 'black'),
                  legend.title = element_text(size = 14,
                                              colour = 'black'),
                  legend.text.align = 0.5,
                  plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black"))

  ggsave(
    filename = sprintf("%s/%s_mito_cutoff.png", output_folder, sample_ID),
    sample_ID_preQC_mito,
    width = 210,
    height = 150,
    units = "mm",
    dpi = 900
  )

}

##############################################################################
# Cell Multiplet plots - scattergram and box plot
##############################################################################
doublet_plots <- function(
  seurat_object,
  sample_ID,
  output_folder,
  ensg_gname_path,
  gene_lower,
  gene_higher_method,
  gene_higher,
  mitochondria_percent,
  nMADs,
  nSD
) {

  # Read in and check the input seurat object
  sample_ID_srt <- read_check_srt(seurat_object, sample_ID, ensg_gname_path)
  # Determine QC parameters
  if (gene_higher_method != "Manual") {
    gene_higher <- out_detect(
      sample_ID_srt,
      sample_ID,
      gene_lower,
      mitochondria_percent,
      outlier_method = gene_higher_method,
      nMADs,
      nSD
    )
  }
  else {
          gene_higher <- as.double(gene_higher)
  }

  # Extract plotting information and create data frame
  sample_ID_preQC_df <- data.frame(sample_ID_srt@meta.data$nFeature_RNA,
                               sample_ID_srt@meta.data$nCount_RNA,
                               sample_ID_srt@meta.data$mito.percent)

  # Rename the columns in the preQC dataframe
  colnames(sample_ID_preQC_df) <- c("gene_count", "UMI_count", "mito_percent")

  # Create new dataframe where all "damaged" cells are removed if thresholds are given
  if (gene_lower != "None"){
    gene_lower <- as.double(gene_lower)
    sample_ID_preQC_df <- sample_ID_preQC_df %>%
      filter(gene_count>gene_lower)
  }
  if (mitochondria_percent != "None"){
    mitochondria_percent <- as.double(mitochondria_percent)
    sample_ID_preQC_df <- sample_ID_preQC_df %>%
      filter(mito_percent<mitochondria_percent)
  }
  sample_ID_nodamaged_df <- sample_ID_preQC_df

  # Add a column to the new dataframe that classifies cells as outliers
  sample_ID_nodamaged_df$outlier_class <- ifelse(sample_ID_nodamaged_df$gene_count>gene_higher,
                                                "Outlier",
                                                "Not outlier")

  # Create the scatter plot flagging cells classified as doublets
  sample_ID_doublet_scatter <- ggplot(sample_ID_nodamaged_df, aes(x = gene_count,
                                                                  y = UMI_count,
                                                                  colour = outlier_class)) +
                               geom_point() +
                               scale_colour_manual(values = c("blue", "#00ff7f"),
                                                    labels = c("Singlet", "Doublet")) +
                               labs(x = "Total number of unique genes detected",
                                    y = "Total number of UMIs detected",
                                    colour = "") +
                               theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                                                         vjust = -5,
                                                                         size = 16),
                                                    axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                                                                vjust = 5,
                                                                                size = 16),
                                                    axis.text.x = element_text(size = 12,
                                                                               colour = 'black'),
                                                    axis.text.y = element_text(size = 12,
                                                                               colour = 'black'),
                                                    legend.text = element_text(size = 12,
                                                                               colour = 'black'),
                                                    legend.title = element_text(size = 14,
                                                                                colour = 'black'),
                                                    legend.text.align = 0.5,
                                                    plot.margin = unit(c(1,1,1,1), "cm")) +
                                theme(panel.background = element_blank(),
                                              axis.line = element_line(colour = "black"))

  ggsave(
    filename = sprintf("%s/%s_doublet_scatter.png",
                        output_folder,
                        sample_ID),
    sample_ID_doublet_scatter,
    width = 210,
    height = 150,
    units = "mm",
    dpi = 900
  )

  # Add in a column to the nodamaged_df to mention the sample name
  sample_ID_nodamaged_df$sample_ID <- sample_ID
  # Create the box plot showing the classified multiplets - for one sample
  sample_ID_doublet_box <- ggplot(sample_ID_nodamaged_df, aes(x = sample_ID,
                                                              y = gene_count,
                                                              fill = "")) +
                           geom_boxplot(outlier.alpha = 0) +
                           scale_fill_manual(values = "#007fff") +
                           geom_point(data = filter(sample_ID_nodamaged_df,
                                                    outlier_class == "Outlier"),
                                      colour = "black",
                                      position = "jitter") +
                           facet_wrap(~ sample_ID, scales = "free_x") +
                           theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 panel.border = element_rect(fill = NA, colour = "black"),
                                 panel.background = element_blank(),
                                 strip.background = element_rect(fill = NA, colour = "black"),
                                 strip.text = element_text(size = 16)) +
                           labs(y = "Total number of unique genes detected",
                                x = "",
                                fill = "") +
                           guides(fill = FALSE,
                                  point = TRUE) +
                           theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                                              vjust = -5,
                                                              size = 16),
                                 axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                                              vjust = 5,
                                                              size = 16),
                                 axis.text.x = element_text(size = 12,
                                                            colour = 'black'),
                                 axis.text.y = element_text(size = 12,
                                                            colour = 'black'),
                                 legend.text = element_text(size = 12,
                                                            colour = 'black'),
                                 legend.title = element_text(size = 14,
                                                            colour = 'black'),
                                 legend.text.align = 0.5,
                                 plot.margin = unit(c(1,1,1,1), "cm"))

  ggsave(
    filename = sprintf("%s/%s_doublet_boxplot.png",
                        output_folder,
                        sample_ID),
    sample_ID_doublet_box,
    width = 210,
    height = 150,
    units = "mm",
    dpi = 900
  )
}

##############################################################################
# Generate the combined QC plot containing all QC metrics with QC parameters
##############################################################################
combined_qc_param <- function(
  seurat_object,
  sample_ID,
  output_folder,
  ensg_gname_path,
  gene_lower,
  gene_higher_method,
  gene_higher,
  mitochondria_percent,
  nMADs,
  nSD
) {

  # Read in and check the input seurat object
  sample_ID_srt <- read_check_srt(seurat_object, sample_ID, ensg_gname_path)
  # Determine QC parameters
  if (gene_higher_method != "Manual") {
    gene_higher <- out_detect(
      sample_ID_srt,
      sample_ID,
      gene_lower,
      mitochondria_percent,
      outlier_method = gene_higher_method,
      nMADs,
      nSD
    )
  }
  # Use the meta.data dataframe from the Seurat object to make the combined QC plot with parameters
  combined_param_plot <- ggplot(sample_ID_srt@meta.data,
                              aes(x = nFeature_RNA,
                                  y = nCount_RNA,
                                  colour = mito.percent)) +
    geom_point() +
    scale_color_gradient(low = 'blue', high = 'red', guide = 'colourbar') +
    labs(x = 'Number of unique genes detected',
         y = 'Number of UMIs detected',
         colour = 'Mitochondrial genes (%)') +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 20),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 20),
          axis.text.x = element_text(size = 14,
                                     colour = 'black'),
          axis.text.y = element_text(size = 14,
                                     colour = 'black'),
          legend.text = element_text(size = 12,
                                     colour = 'black'),
          legend.title = element_text(size = 14,
                                      colour = 'black',
                                      angle = 90,
                                      margin = margin(t = -0.9, unit = 'cm')),
          legend.text.align = 0.5,
          axis.line = element_line(colour = 'black'),
          strip.background = element_rect(fill = 'white', colour = 'black'),
          strip.text = element_text(size = 14)
    ) +
    guides(colour = guide_colorbar(title.position = 'right'))
    # Add layers to the combined QC graph to include the parameters
    if (gene_lower != "None") {
      combined_param_plot <- combined_param_plot +
        geom_vline(xintercept = as.numeric(gene_lower), colour = "limegreen")
    }
    if (gene_higher != "None") {
      combined_param_plot <- combined_param_plot +
        geom_vline(xintercept = as.numeric(gene_higher), colour = "limegreen")
    }
    # Save the plot
    ggsave(sprintf('%s/%s_all_metrics_param.png',
                   output_folder,
                   sample_ID),
           combined_param_plot,
           width = 15,
           height = 10,
           dpi = 600)
}
############################################
# Subset the seurat object based on QC parameters
############################################
sub_srt <- function(
  seurat_object,
  sample_ID,
  output_folder,
  ensg_gname_path,
  gene_lower,
  gene_higher_method,
  gene_higher,
  mitochondria_percent,
  nMADs,
  nSD,
  output_matrix
) {

  # Read in and check the input seurat object
  sample_ID_srt <- read_check_srt(seurat_object, sample_ID, ensg_gname_path)
  # Determine QC parameters
  if (gene_higher != "None"){
    if (gene_higher_method != "Manual") {
      gene_higher <- out_detect(
        sample_ID_srt,
        sample_ID,
        gene_lower,
        mitochondria_percent,
        outlier_method = gene_higher_method,
        nMADs,
        nSD
      )
    }
  }
  # Record the number of cells in the original (input) dataset
  starting_cells = nrow(sample_ID_srt@meta.data)
  cells_remaining = starting_cells

  starting_cells_print_statement = sprintf(
    "%s cells in the input dataset",
    starting_cells
  )
  print(starting_cells_print_statement)

  print("Subsetting the Seurat Object")
  # Subset on lower gene threshold if it is not set to "None"
  if (gene_lower != "None") {
    # Subset the Seurat object
    sample_ID_srt <- subset(sample_ID_srt, subset = nFeature_RNA > as.numeric(gene_lower))
    # Keep record of how many cells were removed and how many cells are remaining
    cells_after_lower_cutoff = nrow(sample_ID_srt@meta.data)
    print_statement_cells_after_lower_cutoff = sprintf("%s cells remaining after applying lower threshold on gene count",
                                                        cells_after_lower_cutoff)
    number_cells_below_threshold = starting_cells - cells_after_lower_cutoff
    print_statement_number_cells_below_threshold = sprintf("%s cells with a gene count below %s",
                                                            number_cells_below_threshold,
                                                            gene_lower)
    percentage_cells_below_threshold = ((starting_cells - cells_after_lower_cutoff)/starting_cells)*100
    print_statement_percentage_cells_below_threshold = sprintf("%s percent of cells with a gene count below %s",
                                                                percentage_cells_below_threshold,
                                                                gene_lower)
    cells_remaining = cells_after_lower_cutoff

    # Use print statements to record the above information in the log file
    print(print_statement_cells_after_lower_cutoff)
    print(print_statement_number_cells_below_threshold)
    print(print_statement_percentage_cells_below_threshold)
  }
  # Subset on mitochondria percentage if it is not set to "None"
  if (mitochondria_percent != "None") {
    # Subset the seurat object
    sample_ID_srt <- subset(sample_ID_srt, subset = mito.percent < as.numeric(mitochondria_percent))

    # Keep record of how many cells were removed and how many cells are remaining
    cells_after_mito_cutoff = nrow(sample_ID_srt@meta.data)
    print_statement_cells_after_mito_cutoff = sprintf("%s cells remaining after applying mitochondria percentage threshold",
                                                      cells_after_mito_cutoff)

    number_cells_above_threshold = cells_remaining - cells_after_mito_cutoff
    print_statement_cells_above_threshold = sprintf("%s cells with a mitochondrial gene percentage higher than %s",
                                                    number_cells_above_threshold,
                                                    mitochondria_percent)
    percentage_cells_above_threshold = ((cells_remaining - cells_after_mito_cutoff)/cells_remaining)*100
    print_statement_percentage_cells_above_threshold = sprintf("%s percent of cells with a mitochondrial percentage above %s",
                                                                percentage_cells_above_threshold,
                                                                mitochondria_percent)

    cells_remaining = cells_after_mito_cutoff

    # Use print statements to record the above information in the log file
    print(print_statement_cells_after_mito_cutoff)
    print(print_statement_cells_above_threshold)
    print(print_statement_percentage_cells_above_threshold)
  }
  # Subset on higher outlier gene count threshold if it is not set to "None"
  if (gene_higher != "None") {
    sample_ID_srt <- subset(sample_ID_srt, subset = nFeature_RNA < as.numeric(gene_higher))

    # Keep record of how many cells were removed and how many cells are remaining
    cells_after_outlier_cutoff = nrow(sample_ID_srt@meta.data)
    print_statement_cells_after_outlier_cutoff = sprintf("%s cells remaining after applying gene count outlier upper threshold",
                                                      cells_after_outlier_cutoff)

    number_cells_above_outlier_threshold = cells_remaining - cells_after_outlier_cutoff
    print_statement_cells_above_outlier_threshold = sprintf("%s cells are detected as outliers, with a gene count above %s",
                                                    number_cells_above_outlier_threshold,
                                                    gene_higher)
    percentage_cells_above_outlier_threshold = ((cells_remaining - cells_after_outlier_cutoff)/cells_remaining)*100
    print_statement_percentage_cells_above_outlier_threshold = sprintf("%s percent of cells are detected as outliers",
                                                                percentage_cells_above_outlier_threshold)

    cells_remaining = cells_after_outlier_cutoff

    # Use print statements to record the above information in the log file
    print(print_statement_cells_after_outlier_cutoff)
    print(print_statement_cells_above_outlier_threshold)
    print(print_statement_percentage_cells_above_outlier_threshold)
  }
  # Save the subset Seurat object based on the default/specified QC parameters
  final_statement = sprintf("Saving subset Seurat object with %s cells to a .rds file",
                            cells_remaining)
  print(final_statement)
  saveRDS(sample_ID_srt,
          file = sprintf('%s/%s_subset_seurat.rds', output_folder, sample_ID))

  # Save the counts in the subset seurat object as an mtx count matrix if output_matrix == TRUE
  if (output_matrix == TRUE) {
    srt_to_mtx(sample_ID_srt, sample_ID, output_folder)
  }
}

####################################################################
# R Function to extract cellQC metrics - prior to subsetting
####################################################################
extract_metrics <- function(
  seurat_object,
  sample_ID,
  output_folder,
  ensg_gname_path,
  gene_lower,
  gene_higher_method,
  gene_higher,
  mitochondria_percent,
  nMADs,
  nSD
){
  # Read in the seurat object
  sample_ID_srt <- read_check_srt(seurat_object, sample_ID, ensg_gname_path)
  # Determine QC parameters
  if (gene_higher_method != "Manual") {
    gene_higher <- out_detect(
      sample_ID_srt,
      sample_ID,
      gene_lower,
      mitochondria_percent,
      outlier_method = gene_higher_method,
      nMADs,
      nSD
    )
  }
  print("Extracting cell quality control metrics")
  # Extract plotting information and create data frame
  sample_ID_preQC_df <- data.frame(sample_ID_srt@meta.data$nFeature_RNA,
                               sample_ID_srt@meta.data$nCount_RNA,
                               sample_ID_srt@meta.data$mito.percent)

  # Rename the columns in the preQC dataframe
  colnames(sample_ID_preQC_df) <- c("gene_count",
                                    "UMI_count",
                                    "mito_percent")

  # Use lower gene count threshold to classify cells as damaged
  if (gene_lower != "None"){
    gene_lower <- as.numeric(gene_lower)
    sample_ID_preQC_df$lower_threshold <- ifelse(sample_ID_preQC_df$gene_count<gene_lower,
                                                  "Damaged",
                                                  "Healthy")
  }
  # Use mitochondrial percentage threshold to classify cells as damaged
  if (mitochondria_percent != "None"){
    mitochondria_percent <- as.numeric(mitochondria_percent)
    sample_ID_preQC_df$mitochondria_threshold <- ifelse(sample_ID_preQC_df$mito_percent>mitochondria_percent,
                                                        "Damaged",
                                                        "Healthy")
  }
  # Add a column to the new dataframe that classifies cells as outliers
  if (gene_higher != "None"){
    gene_higher <- as.numeric(gene_higher)
    sample_ID_preQC_df$outlier_class <- ifelse(sample_ID_preQC_df$gene_count>gene_higher,
                                                  "Multiplet",
                                                  "Singlet")
  }

  # Write the dataframe to a tab delimited file
  print("Writing cell quality control metrics to dataframe")
  write_tsv(
    sample_ID_preQC_df,
    path = sprintf("%s/%s_cell_metrics_dataframe.tsv",
                    output_folder,
                    sample_ID)
  )
}

####################################################################
# R Function to generate all cellQC default plots that are included
####################################################################
default_plots <- function(
  seurat_object,
  sample_ID,
  output_folder,
  ensg_gname_path,
  gene_lower,
  gene_higher_method,
  gene_higher,
  mitochondria_percent,
  nMADs,
  nSD
){
  print("Creating default plots for cell quality control")
  # Generate violin plots
  vln_srt(
    seurat_object,
    sample_ID,
    output_folder,
    ensg_gname_path
  )
  # Generate the combined plots
  combined_qc(
    seurat_object,
    sample_ID,
    output_folder,
    ensg_gname_path
  )
  # Generate all plots representing damaged cells
  if (gene_lower != "None" & mitochondria_percent != "None"){
    damaged_scatter(
      seurat_object,
      sample_ID,
      output_folder,
      ensg_gname_path,
      gene_lower,
      mitochondria_percent
    )
  }
  # Generate all plots representing cell multiplets
  if (gene_higher != "None"){
    doublet_plots(
      seurat_object,
      sample_ID,
      output_folder,
      ensg_gname_path,
      gene_lower,
      gene_higher_method,
      gene_higher,
      mitochondria_percent,
      nMADs,
      nSD
    )
  }
}


########################################################
# R function to run QC on single-cell data using Seurat
########################################################
seurat_qc <- function(
  file,
  sample_ID,
  output_folder,
  generate_seurat_object,
  subset_seurat_object,
  generate_default_plots,
  gene_column,
  input_seurat_object,
  ensg_gname_path,
  gene_lower,
  gene_higher_method,
  gene_higher,
  mitochondria_percent,
  nMADs,
  nSD,
  extract_cell_metrics,
  output_matrix
){
  if (generate_seurat_object == TRUE) {
    seurat_object <- gen_srt(
      file,
      sample_ID,
      output_folder,
      ensg_gname_path,
      gene_column
    )
  } else {
    seurat_object <- file
  }
  if (subset_seurat_object == TRUE) {
    sub_srt(
      seurat_object,
      sample_ID,
      output_folder,
      ensg_gname_path,
      gene_lower,
      gene_higher_method,
      gene_higher,
      mitochondria_percent,
      nMADs,
      nSD,
      output_matrix
    )
  }
  if (generate_default_plots == TRUE){
    default_plots(
      seurat_object,
      sample_ID,
      output_folder,
      ensg_gname_path,
      gene_lower,
      gene_higher_method,
      gene_higher,
      mitochondria_percent,
      nMADs,
      nSD
    )
  }
  # Extract the cell metrics into a tab separated file if required
  if (extract_cell_metrics == TRUE){
    extract_metrics(
      seurat_object,
      sample_ID,
      output_folder,
      ensg_gname_path,
      gene_lower,
      gene_higher_method,
      gene_higher,
      mitochondria_percent,
      nMADs,
      nSD
    )
  }
}

args = commandArgs(trailingOnly = TRUE)
seurat_qc(
  args[1],
  args[2],
  args[3],
  args[4],
  args[5],
  args[6],
  args[7],
  args[8],
  args[9],
  args[10],
  args[11],
  args[12],
  args[13],
  args[14],
  args[15],
  args[16],
  args[17]
)
