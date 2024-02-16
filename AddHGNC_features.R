

# required packages
library(tidyverse)

# Function to add the HGNC symbols to the feature file with the t2g as reference
AddHGNC <- function(
    feature_file, # Path to the features file with ENSG feature names
    t2g, # Path to the transcripts to genes mapping file to be used to add HGNC to the ENSG features file
    adjusted_features_file_output, # Where to store the adjusted features file
    adjusted_t2g_output # Where to store the adjusted t2g file (adjusted if t2g file has missing values)
){ # 4 parameters and no defaults
  # Step one: Read in the files
  FeatureFromBustools <- read_tsv(
    feature_file,
    col_names = FALSE # There are no column names so this parameter is essential to specify
  ) # Remember this will have to be given by the user
  # The transcripts to genes mapping file
  T2GFile <- read_tsv(
    t2g,
    col_names = FALSE # There are no column names so this parameter is essential to specify
  ) # The path to this file also has to be supplied by the user
  # Step two: Create a new column in the t2g file to ensure there are no missing HGNC symbols
  T2GFile <- mutate(
    T2GFile,
    X4 = case_when( # We have to create a new column otherwise specifiying the same columns overwrites the current data
      is.na(X3) ~ X2, # This doesn't really seem like it would work
      !is.na(X3) ~ X3
    )
  ) # new file now contains four columns 
  # Now I just need to delete the X3 column (Because X4 is now the new X3)
  # If there was any missing data X4 will solve the problem
  # If there wasn't any missing data X4 will be a duplicate of X3
  T2GFile <- dplyr::select(
    T2GFile,
    c(
      X1,
      X2,
      X4
    )
  ) # now this file will have 3 columns (ENST, ENSG, HGNC)
  # Step three: Save the adjusted t2g file (Or it would just be a duplicate)
  t2g_output_filename <- sprintf(
    "%s/adjusted_t2g.txt",
    adjusted_t2g_output
  )
  write_tsv(
    x = T2GFile,
    file = t2g_output_filename
  )
  # Step four: select X2 and X4 from t2g and remove duplicates
  # Select the gene columns we want to work with
  T2GFile <- dplyr::select(
    T2GFile,
    c(
      X2,
      X4
    )
  )
  # Rename the column names so that it is easier to work with
  colnames(T2GFile) <- c("ENSG", "HGNC")
  # We should also rename the FeatureFromBustools column so that it matches
  colnames(FeatureFromBustools) <- c("ENSG")
  # Now we need to remove duplicates from the T2GFile dataframe
  T2GFile <- dplyr::distinct(
    T2GFile
  )
  # Step five: add the HGNC symbols to the BUStools feature ENSG feature file
  FeatureFromBustoolsHGNC <- left_join(
    x = FeatureFromBustools,
    y = T2GFile,
    by = "ENSG"
  )
  # Step six: Save the new features file 
  adjusted_features_file_name <- sprintf(
    "%s/features.tsv",
    adjusted_features_file_output
  )
  # Save it as a tsv and use write_tsv to do so
  write_tsv(
    x = FeatureFromBustoolsHGNC,
    file = adjusted_features_file_name, # This path will need to be specified by the user as well
    col_names = FALSE
  )
  # And that's it - the rest happens in bash
}

# Now allow this function to be run through bash
args = commandArgs(trailingOnly = TRUE)
AddHGNC(
  args[1],
  args[2],
  args[3],
  args[4]
)