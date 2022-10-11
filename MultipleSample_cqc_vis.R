
# required packages
library(tidyverse)
library(Seurat)
library(patchwork)
library(gridExtra)
library(Matrix)
library(stats)

############################# Function for editing plots ##########
# Function for editing plots
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 14),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 14),
          axis.text.x = element_text(size = 12,
                                     colour = 'black'),
          axis.text.y = element_text(size = 12,
                                     colour = 'black'),
          legend.text = element_text(size = 12,
                                     colour = 'black'),
          legend.title = element_text(size = 14,
                                      colour = 'black'),
          strip.text = element_text(
            size = 14,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

####################################################################################
# R function to create visualisations to view QC metrics across multiple samples
####################################################################################

many_samples_vis <- function(
        all_input_directories, # As string since it is coming from bash
        all_sample_IDs, # As string since it is coming from bash
        main_output_directory, # Where to store the combined sample visualisations
        output_prefix # What the prefix of the output visualisations should be
) {
  # First thing we need to do is convert the input_directory string to a vector
  all_input_directories_vector <- strsplit(
    all_input_directories,
    "," # Will be combined with , in the outer python function
  )[[1]] # The first element of the list is where out vector will be stored
  # Convert the samples_IDs string to a vector
  all_sample_IDs_vector <- strsplit(
    all_sample_IDs,
    ","
  )[[1]]
  # Read in the tables for each of the samples (will be stored as a csv) - will need to use a for loop
  cellQCsDfList <- list()
  for ( i in seq(1, length(all_sample_IDs_vector)) ) { # Best to always use indices
    # What is the name used to store the CellCountInfo Df
    DfName <- sprintf(
      "%s/%s_CellCountInfo.csv",
      all_input_directories_vector[i],
      all_sample_IDs_vector[i]
    )
    # Read in the Df using read_csv (since it was saved as csv file)
    CellInfoDf <- read_csv(
      DfName
    )
    # Store that Df in the list opened outside the loop
    cellQCsDfList <- c(
      cellQCsDfList,
      list(CellInfoDf)
    )
  }
  # Combine all of those DFs stored in the list into one major Df
  cellQCsCombinedDf <- bind_rows(
    cellQCsDfList
  )
  # Now wrangle the dataframe so the proportions can be calculated and it can be plotted
  cellQCsGatheredDf <- pivot_longer(
    cellQCsCombinedDf,
    cols = c(
      BelowGeneCount,
      AboveMt,
      UpperOutlier,
      CellsLeft
    ),
    names_to = "CellClass",
    values_to = "CellNumber"
  )
  # Calculate the proportions using mutate
  cellQCsGatheredDf <- mutate(
    cellQCsGatheredDf,
    CellProps = (CellNumber/Starting)*100
  )
  # Now we can plot
  cellQCsPropsStackBar <- ggplot(
    cellQCsGatheredDf,
    aes(
      x = Sample,
      y = CellProps,
      fill = CellClass
    )
  ) +
    geom_bar(
      stat = "identity",
      width = 0.7,
      colour = "black"
    ) +
    labs(
      x = "",
      y = "Proportion of Cells (%)",
      fill = "Category"
    ) +
    scale_fill_manual(
      values = c(
        "#ff007f",
        "#fff000",
        "blue",
        "#7fff00"
      ),
      labels = c(
        "Mitocondrial Gene\nPercentage > 10%",
        "Total Gene Count < 200",
        "High quality Cells",
        "Potential multiplets"
      )
    )
  cellQCsPropsStackBar <- edit_plots(cellQCsPropsStackBar)
  # Now save the plot
  PlotName <- sprintf(
    "%s/%s_CombinedCellQCVis.png",
    main_output_directory,
    output_prefix
  )
  # Use ggsave
  ggsave(
    filename = PlotName,
    plot = cellQCsPropsStackBar,
    width = 8,
    height = 9,
    dpi = 900
  )
  # Add a return function just so we can assess the output
  # Since it doesn't look as expected - return the uncombined list - mistake was a dumb one
  #return(cellQCsPropsStackBar) # This works, now we just need to adjust labels and save
}

args = commandArgs(trailingOnly = TRUE)
many_samples_vis(
  args[1],
  args[2],
  args[3],
  args[4]
)
