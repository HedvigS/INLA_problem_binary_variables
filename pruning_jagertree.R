#This script counts how often words occur in language names in Glottolog. It only takes into account the "main name" and does not include dialects, proto-languages or other non-languages (artifical, unnatested etc).

#installing and loading packages
if (!suppressPackageStartupMessages(require("pacman"))) { install.packages("pacman") } #if pacman isn't already installed, install it.

pacman::p_load(
tidyverse, 
ape, 
phytools
  )

options(tidyverse.quiet = TRUE)

#script written by Sam Passmore and Hedvig Skirgård

cat("Pruning the global tree from Jäger to the dataset.\n")

#### Inputs ####
#Glottolog-cldf table to have a match for all dialects to their language parent. Note that the particular dialects may differ from the dialects in GB which is why we cann't use the language table from the the dataset-cldf relase
glottolog_fn <- "data/glottolog_df.tsv"
if (!file.exists(glottolog_fn)) {
source("make_glottolog-cldf_table.R") }
glottolog_df <- read_tsv(glottolog_fn, col_types = cols()) %>% 
  dplyr::select(Longitude, Latitude, Language_ID, Language_level_ID) %>% 
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID))

taxa_pairing <- read.csv('data/taxa.csv') %>% 
  dplyr::rename(Language_ID = glottocode) %>% 
  left_join(glottolog_df, by = "Language_ID")

jaeger_tree <- read.tree('data/world.tre')

languages <- read_tsv("data/Sahul_structure_wide.tsv",col_types = cols()) %>% 
  dplyr::select(Language_ID = glottocode) #this column is already aggregated for dialects in make_wide.R

#### Tree - Data pairs ####
coverage = sum(languages$Language_ID %in% taxa_pairing$Language_ID) / nrow(languages) 
cat("Tips in the global Jaeger tree  can be matched to", round(coverage, 2) * 100, "% of the langauges in the dataset.\n")

#### Subset to the dataset langauges ####
in_tree <- jaeger_tree$tip.label %>% 
  as.data.frame() %>% 
  dplyr::rename(taxon_full = ".") %>% #making a data frame with a column representing all the tip labels, in the right order
  mutate(taxon = str_extract(taxon_full, "[^.]+$")) %>% #extracting the part of the tip label that matches to the taxa file
  left_join(taxa_pairing, by = "taxon") %>% #pairing with the taxa file
  dplyr::select(taxon_full, Language_ID = Language_level_ID) %>% #only keeping the necessary cols for pruning and matching with GB
  inner_join(languages, by = "Language_ID") %>% #prune to only tips which are also in GB
  group_by(Language_ID) %>%
  sample_n(size = 1)  #if there is more than one tip with the same language ID (duplicates or dialects), randomly remove all but one

jaeger_pruned = keep.tip(jaeger_tree, in_tree$taxon_full)

#### Relabel taxa in tree ####
#taxons are matched to the glottocode of themselves or their parent which has the level "language" in glottolog

jaeger_pruned$tip.label <- jaeger_pruned$tip.label %>% 
  as.data.frame() %>% 
  dplyr::rename(taxon = ".") %>% 
  mutate(taxon = str_extract(taxon, "[^.]+$")) %>% 
  left_join(taxa_pairing, by = "taxon") %>%
  dplyr::select(Language_level_ID) %>% 
  .[,1]

#scaling tree
scale = 1
jaeger_pruned$edge.length = jaeger_pruned$edge.length/max(nodeHeights(jaeger_pruned)[,2])*scale

write.tree(jaeger_pruned, "data/jaeger_pruned.tree")

cat("Pruned Jager tree created with", length(jaeger_pruned$tip.label), "tips.\n")