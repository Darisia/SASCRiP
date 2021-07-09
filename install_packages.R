# R script containing code to install required packages

# First check if the package is installed - if yes check the version - if the version
# .. is not correct return a warning

if ("tidyverse" %in% rownames(installed.packages())){
  print("tidyverse is already installed")
} else {
  install.packages("tidyverse")
}

if ("Seurat" %in% rownames(installed.packages())){
  print("Seurat is already installed")
} else {
  install.packages("Seurat")
}
