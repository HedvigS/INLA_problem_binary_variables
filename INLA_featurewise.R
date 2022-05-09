#This is a script for running binomial INLA over a set of linguistic features, with phylo and spatial effects.

set.seed(67)

#set this as 1 if you're just running this script on 50 lgs over 3 features to debug. Otherwise set to 0.
debug_run = 0

#installing and loading packages
source("requirements.R")

#for detect_coderbias.R and spatiophylogenetic_jaegermodel.R
kappa = 2 # smoothness parameter as recommended by Dinnage et al. (2020)
sigma = c(1, 1.15) # Sigma parameter. First value is not used. 

# load variational covariance matrix function taken from geoR::varcov_spatial
source("varcov_spatial.R")

# Check that INLA is installed
source("install_inla.R")

#dirs
if (!dir.exists("output" )) {
  dir.create(  "output" )
}

if(debug_run == 1){
  OUTPUTDIR  <- file.path("output", "results_debug/")
} else{
  OUTPUTDIR <- file.path("output", "results/")
}

if (!dir.exists(  OUTPUTDIR )) {
  dir.create(  OUTPUTDIR )
}

if (!dir.exists(file.path(  OUTPUTDIR , "phylo_only/"))) {
  dir.create(file.path(  OUTPUTDIR , "phylo_only/"))
  dir.create(file.path(  OUTPUTDIR , "spatial_only/"))
  dir.create(file.path(  OUTPUTDIR , "autotyp_area_only/"))
  dir.create(file.path(  OUTPUTDIR , "dual_process_rdata/"))
  dir.create(file.path(  OUTPUTDIR , "trial_process_rdata/"))
  }		

sink(file = file.path(  OUTPUTDIR , "INLA_featurewise_log.txt"), split = T)

cat("Starting INLA runs at", as.character(Sys.time()), ".\n")

#### Functions ####

cov2precision = function(spatial_covar_mat){
  spatial_covar_mat = spatial_covar_mat / exp(determinant(spatial_covar_mat)$modulus[1] /
                                                nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  spatial_prec_mat
}

# useful objects
join_columns = c("2.5%", "50%", "97.5%", "Feature_ID", 
                 "effect", "waic", "model")


#### Main Analyses ####

cat("#### Building Jaeger tree models ####\n")

#reading in data
Sahul_data_fn <- "data/Sahul_structure_wide_imputed.tsv"
if (!file.exists(Sahul_data_fn)) {
  source("Impute_missing_values.R")
}
data <- read_tsv(Sahul_data_fn ,col_types = cols()) %>% 
  dplyr::rename(Language_ID = ID) #this column is already aggregated for dialects in make_wide.R

#subset GB to test code for debugging
if(debug_run == 1){
data <- data[1:150,]
}

#### Inputs ####
# language metadata
if (!file.exists("data/glottolog_AUTOTYP_areas.tsv")) { source("assigning_AUTOTYP_areas.R") }		
autotyp_area <- read.delim("data/glottolog_AUTOTYP_areas.tsv", sep = "\t") %>%
  dplyr::select(Language_ID, AUTOTYP_area)

glottolog_df_fn = "data/glottolog_df.tsv"
if (!file.exists(glottolog_df_fn)) { source("make_glottolog-cldf_table.R") }		

languages <- read.delim(glottolog_df_fn, sep = "\t") %>%		
  dplyr::select(Language_ID, Family_ID, Name, Longitude, Latitude, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  inner_join(dplyr::select(data, "Language_ID"), by = "Language_ID") %>% 
  mutate(Longitude = round(Longitude, 3)) %>% # let's cut down the precision of the lat/long to make the process go quicker. See stack exchange thread where they say "The third decimal place is worth up to 110 m: it can identify a large agricultural field or institutional campus." https://gis.stackexchange.com/questions/8650/measuring-accuracy-of-latitude-and-longitude
  mutate(Latitude = round(Latitude, 3)) %>% 
  left_join(autotyp_area, by = "Language_ID")

# trees
tree_filename = 'data/jaeger_pruned.tree'
if (!file.exists(tree_filename)) { source("pruning_jagertree.R") }		
phylogenetic_tree = read.tree(tree_filename)

# Subset GB and languages to Jaeger set
data  <- data[data$Language_ID %in% phylogenetic_tree$tip.label,]
languages = languages[languages$Language_ID %in% data $Language_ID,]

# prune tree to dataset
taxa = data$Language_ID
phylogenetic_tree = keep.tip(phylogenetic_tree, tip = taxa)

#### Parameters ####
#### Spatial Jittering ####
## There are some number of languages that have identical spatial coordinates, which we cannot allow for the spatial analysis.
## I have jittered the coordinates that are identical
## Jittering moves locations randomly by about 0.3 and 1 degree in Longitude & Latitude
duplicate_coords = languages[duplicated(languages[,c("Longitude", "Latitude")]) | 
                               duplicated(languages[,c("Longitude", "Latitude")], 
                                          fromLast = TRUE),"Language_ID"]
duplicate_rowid = languages$Language_ID %in% duplicate_coords
languages$Latitude[duplicate_rowid] = jitter(languages$Latitude[duplicate_rowid], 
                                             factor = 1)
languages$Longitude[duplicate_rowid] = jitter(languages$Longitude[duplicate_rowid], 
                                              factor = 1)

#### Phylogenetic covariance matrix ####
cat("Calculating the phylogenetic variance covariance matrix.\n")

phylo_covar_mat <- ape::vcv(phylogenetic_tree)
phylo_covar_mat <- phylo_covar_mat / max(phylo_covar_mat)
# The diagonal of phylo_covar_mat should inform our prior
phylo_prec_mat = cov2precision(phylo_covar_mat)
#### Spatial covariance matrix ####

cat("Calculating the spatial variance covariance matrix.\n")
## Ensure the order of languages matches the order within the phylogeny
languages = languages[order(match(languages$Language_ID, rownames(phylo_prec_mat))),]

spatial_covar_mat = varcov.spatial(languages[,c("Longitude", "Latitude")], 
                                   cov.pars = sigma, kappa = kappa)$varcov
dimnames(spatial_covar_mat) = list(languages$Language_ID, languages$Language_ID)

spatial_prec_mat = cov2precision(spatial_covar_mat)


#### Set up model priors ####

## Taken from Dinnage et al. (2020):
### As in the manuscript, we use a “PC” prior, which stand for “Penalizing Complexity”.
### This is a standard prior developed by the developers of INLA, which is “weakly informative”.
### It places slightly more of the prior probability density on values close to zero,
### but has a “long tail”, which allows the data to push the parameter away from zero if
### there is good evidence (i.e. the likelihood of the data is higher).
### The PC prior has two parameters, p1 and p2: p2 is the proportion of the prior probability
### density that falls above values greater than p1.
###
### Given we are modelling (somewhat) Gaussian data that we have standardised to a variance of 1,
### we would not expect any random factor to have a variance greater than one.
### So we will set our prior to only have about 10% of its prior probability density
### above 1. We will guestimate the total variance possibly explained by the
### phylogenetic effect based on it’s diagonal entries, which are all equal
### to 2.7373648. So to get to 1, the scaling factor would have to be about 0.36

pcprior = list(prec =list(prior="pc.prec", param = c(1, 0.1)))

kappa = 2 # smoothness parameter as recommended by Dinnage et al. (2020)
sigma = c(1, 1.15) # Sigma parameter. First value is not used. 


## Adding random effect ids
df = data %>%
  left_join(tibble(Language_ID = rownames(phylo_prec_mat),
                  phy_id_generic = 1:nrow(phylo_prec_mat),
                  phy_id_iid_model = 1:nrow(phylo_prec_mat),
                  spatial_id_generic = 1:nrow(spatial_prec_mat),
                  spatial_id_iid_model = 1:nrow(spatial_prec_mat)), 
            by = "Language_ID") %>% 
  left_join(languages,  by = "Language_ID") %>% 
  rename(AUTOTYP_area_id_iid_model = AUTOTYP_area)

#################
###INLA LOOPS####
#################

#features to loop over
features <- data %>% 
  dplyr::select(-Language_ID) %>% 
  colnames() 

#subsetting for debugging code swiftly
if(debug_run == 1) {
features <- features[3:7]
}

cat("#### Phylogenetic only model ####\n")

#make empty df to bind to
df_phylo_only <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(df_phylo_only) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "model", "waic", "marginals.hyperpar.phy_id_iid_model", "marginals.hyperpar.phy_id_generic") 
df_phylo_only$`2.5%` <- as.numeric(df_phylo_only$`2.5%`)
df_phylo_only$`50%` <- as.numeric(df_phylo_only$`50%`)
df_phylo_only$`97.5%` <- as.numeric(df_phylo_only$`97.5%`)
df_phylo_only$Feature_ID <- as.character(df_phylo_only$Feature_ID)
df_phylo_only$model <- as.character(df_phylo_only$model)
df_phylo_only$effect <- as.character(df_phylo_only$effect)
df_phylo_only$waic <- as.numeric(df_phylo_only$waic)
df_phylo_only$marginals.hyperpar.phy_id_iid_model <- as.list(df_phylo_only$marginals.hyperpar.phy_id_iid_model)
df_phylo_only$marginals.hyperpar.phy_id_generic <- as.list(df_phylo_only$marginals.hyperpar.phy_id_generic)

index <- 0

cat("Starting INLA phylo-only featurewise runs at", as.character(Sys.time()), ".\n")

for(feature in features){
  
  #feature <- features[149]
  
  index <- index + 1   
  cat(paste0("# Running the phylo-only model on feature ", 
             feature, 
             ". This feature has index ", index, ". That means I'm ", 
             round(index/length(features) * 100, 
                   2), 
             "% done.\n"))
  

  if(feature %in% colnames(df)){
  
  cat("Feature is in df, progressing..\n")
  formula <- eval(substitute(this_feature ~
                               f((phy_id_generic), 
                                 model = "generic0",
                                 Cmatrix = phylo_prec_mat,
                                 constr = TRUE, 
                                 hyper = pcprior) + 
                               f(phy_id_iid_model,
                                 model = "iid", 
                                 hyper = pcprior), list(this_feature=as.name(feature))))
  
  output <- try(expr = {
    inla(formula = formula,
                                 control.compute = list(waic=TRUE, dic = FALSE, mlik = FALSE, config = TRUE),
                                 control.inla = list(tolerance = 1e-6, h = 0.001),
                                 control.predictor = list(compute=TRUE, link=1), 
                                 control.family = list(control.link=list(model="logit")),  
                                 data = df,family = "binomial")})

  
  if (class(output) != "try-error") {
  
    cat("I've finished running inla for feature ", feature, ", writing RDS now and extracting effects.\n")
suppressWarnings(  saveRDS(output, file = paste0(OUTPUTDIR, "phylo_only/phylo_only_", feature, ".rdata")) )
#Don't be alarmed by the suppress warnings. saveRDS() is being kind and reminding us that the package stats may not be available when loading. However, this is not a necessary warning for us so we've wrapped saveRDS in suppressWarnings.

#pulling out phy_id_generic effect
#if the hessian has negative eigenvalues, then the hyperpar will contain inf values and the extract won't work, therefore there's an if statement testing for this.

if(!(Inf %in% output$marginals.hyperpar$`Precision for phy_id_generic`[,2])){

phylo_effect_generic = inla.tmarginal(function(x) 1/x,
                                output$marginals.hyperpar$`Precision for phy_id_generic`,
                                method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)

df_phylo_only_generic  <- phylo_effect_generic %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
  mutate(Feature_ID = feature) %>% 
  mutate(effect = "phylo_only_generic") %>% 
  mutate(model = "phylo_only") %>% 
  mutate(waic = output$waic$waic)  %>% 
  mutate(marginals.hyperpar.phy_id_generic = output$marginals.hyperpar[1])
}else{
  df_phylo_only_generic <- tibble(
    "2.5%" = c(NA),
    "50%" =c(NA),
    "97.5%" =c(NA)) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "phylo_only_generic") %>% 
    mutate(model = "phylo_only") %>% 
    mutate(waic = output$waic$waic)  %>% 
  mutate(marginals.hyperpar.phy_id_generic = output$marginals.hyperpar[1])
}

#pulling out phy_id_iid_model effect
#if the hessian has negative eigenvalues, then the hyperpar will contain inf values and the extract won't work, therefore there's an if statement testing for this.

if(!(Inf %in% output$marginals.hyperpar$`Precision for phy_id_iid_model`[,2])){
phylo_effect_iid_model = inla.tmarginal(function(x) 1/x,
                                        output$marginals.hyperpar$`Precision for phy_id_iid_model`,
                                        method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .) 

df_phylo_only_iid_model <- phylo_effect_iid_model %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
  mutate(Feature_ID = feature) %>% 
  mutate(effect = "phylo_only_iid_model") %>% 
  mutate(model = "phylo_only") %>% 
  mutate(waic = output$waic$waic)  %>% 
  mutate(marginals.hyperpar.phy_id_iid_model = output$marginals.hyperpar[2])
} else{
  df_phylo_only_iid_model<- tibble(
    "2.5%" = c(NA),
    "50%" =c(NA),
    "97.5%" =c(NA)) %>%  
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "phylo_only_iid_model") %>% 
    mutate(model = "phylo_only") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.phy_id_iid_model = output$marginals.hyperpar[2])
}

df_phylo_only <- df_phylo_only  %>% 
  full_join(df_phylo_only_iid_model, 
            by = c(join_columns, "marginals.hyperpar.phy_id_iid_model")) %>% 
  full_join(df_phylo_only_generic, 
            by = c(join_columns, "marginals.hyperpar.phy_id_generic"))

  }
  }else{
    cat("Feature", feature, "wasn't in the dataframe. Moving on.")
  }
rm(output)  
}
cat("All done with the phylo only model, 100% done!")

df_phylo_only %>% write_tsv(file = file.path(OUTPUTDIR, "df_phylo_only.tsv"))
df_phylo_only %>% saveRDS(file = file.path(OUTPUTDIR, "df_phylo_only.Rdata"))

sink(file = NULL)