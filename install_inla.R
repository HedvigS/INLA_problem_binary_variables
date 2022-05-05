#installing and loading packages
if (!suppressPackageStartupMessages(require("pacman"))) { install.packages("pacman") } #if pacman isn't already installed, install it.

p_load(rlang)

# Check that INLA is installed
if (!rlang::is_installed("INLA")) { 
  cat("INLA wasn't installed, installing now.\n") 
  # 1. Install BiocManager using: 
  pacman::p_load("BiocManager")
  
  # 2. Install INLA dependencies with BiocManager using: 
  BiocManager::install(c("graph","Rgraphviz",
                         "rgdal",
                         "rgl",
                         "spdep"))
  
  # 3. Update foreach (although it unclear how vital this step was) using: 
  #install.packages("foreach")
  
  # 4. Install INLA using: 
  # NOTE: This is a big download
  install.packages("INLA", repos=c(getOption("repos"), 
                                   INLA="https://inla.r-inla-download.org/R/stable"), 
                   dep=TRUE)
  } else {
  cat("Great, INLA was already installed, loading now.\n") 
}
suppressPackageStartupMessages(
  library(INLA, quietly = T, warn.conflicts = F, verbose = F)
)

inla.setOption(inla.mode="experimental") 