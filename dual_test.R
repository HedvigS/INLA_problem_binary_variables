## Single INLA test
## Binomial INLA test
source("requirements.R")

# I haven't gone through to check which libraries are loaded
# within requirements and had already added these.
# Can clean up later. 
suppressPackageStartupMessages({
  library(spdep)
  library(phytools)
  library(geiger)
  library(ape)
  library(INLA)
  library(caper)
  library(dplyr)
  library(assertthat)
  library(geoR)
})

## Parameters
args = commandArgs(trailingOnly=TRUE)

lambda = args[1]

## functions
cov2precision = function(spatial_covar_mat){
  spatial_covar_mat = spatial_covar_mat / 
    exp(determinant(spatial_covar_mat)$modulus[1] /
          nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  spatial_prec_mat
}

get_lambda_inla = function(fit, effect){
  hyper_summary = fit$summary.hyperpar
  eff_1 = hyper_summary[effect, "mean"]
  eff_2 = sum(hyper_summary[,"mean"])
  binomial_constant = pi^2/3
  (1 / eff_1) / 
    ((1/eff_1) + (1/eff_2) + (1/binomial_constant))
}

# grambank metadata
grambank_metadata = read_csv("data/languages.csv",
                             col_types = cols()) %>%		
  dplyr::select(Language_ID = Language_level_ID, 
                Name, 
                Longitude, 
                Latitude) %>% 
  distinct(Language_ID, .keep_all = T)

# jager tree
tree = read.tree('./data/jaeger_pruned.tree')

## Subset grambank to all of those in Jager tree
keep_languages = grambank_metadata$Language_ID %in% tree$tip.label
grambank_metadata = grambank_metadata[keep_languages,]

# Sort Grambank to match tree tips
grambank_metadata = 
  grambank_metadata[match(tree$tip.labe,
                          grambank_metadata$Language_ID),]


# Use real longitude and latitude from Grambank
longitude = grambank_metadata$Longitude
latitude = grambank_metadata$Latitude

model_data = data.frame(longitude = longitude,
                        latitude = latitude,
                        phy_id_int = 1:nrow(grambank_metadata),
                        spat_id_int =  1:nrow(grambank_metadata),
                        error_int = 1:nrow(grambank_metadata),
                        glottocodes1 = grambank_metadata$Language_ID,
                        glottocodes2 = grambank_metadata$Language_ID,
                        glottocodes3 = grambank_metadata$Language_ID)

## Make matrices
#### Spatial matrix
cat("Calculating the spatial variance covariance matrix.\n")
## Ensure the order of languages matches the order within the phylogeny
spatial_covar_mat = varcov.spatial(model_data[,c("longitude", "latitude")], 
                                   cov.pars = sigma, kappa = kappa)$varcov
dimnames(spatial_covar_mat) = list(model_data$glottocodes, model_data$glottocodes)
spatial_prec_mat = cov2precision(spatial_covar_mat)

## Make Phylogenetic matrix
phylo_covar_mat <- ape::vcv(tree)
phylo_prec_mat = cov2precision(phylo_covar_mat)

# Priors
pcprior_phy = list(prec = list(
  prior="pc.prec",
  param = c(1, 0.1)) # probability that lambda is 0.1 is 10%
)

output_list = list()
iter = 20
for(i in 1:iter){
  print(i)
  y = rTraitDisc(
    geiger::rescale(tree,
            lambda,
            model = "lambda"),
    k = 2,
    freq = 0.5,
    states = 0:1
  )
  
  model_data$y = as.numeric(y) - 1
  
  print("fitDiscrete...")
  pagels_lambda = fitDiscrete(tree, 
                              factor(y), 
                              transform = "lambda")
  
  print("INLA...")
  inla_model = inla(y ~ 
         f(phy_id_int,
           model = "generic0",
           Cmatrix = phylo_prec_mat,
           constr = TRUE,
           hyper = pcprior_phy) + 
         f(spat_id_int,
           model = "generic0",
           Cmatrix = spatial_prec_mat,
           constr = TRUE,
           hyper = pcprior_phy) +
         f(error_int,
           model = "iid",
           constr = TRUE,
           hyper = pcprior_phy) , 
       data = model_data)
  
  print("brms...")
  brms_model <- brm(
    y ~ 1 + 
      (1|gr(glottocodes1, cov = spatial)) + 
      (1|gr(glottocodes2, cov = phylogeny)) + 
      (1|glottocodes3), 
    data = model_data, 
    family = bernoulli(), 
    data2 = list(spatial = spatial_covar_mat,
                 phylogeny = phylo_covar_mat),
    prior = c(
      prior(normal(0, 50), "Intercept"),
      prior(student_t(3, 0, 20), "sd")
    ),
    chains = 1
  )
  
  print("Phylo D...")
  phylo_d_results = phylo.d(data = model_data,
                              names.col = glottocodes, 
                              phy = tree,
                              binvar = y)
  
  output_list[[i]] = list(y = y,
                          pagels_lambda = pagels_lambda,
                          inla_model = lambda_model,
                          brms_model = brms_model,
                          phylo_d = phylo_d_results)
}

saveRDS(output_list, file = 
          paste0(
            "temp_scripts_for_meetings/spatial",
            lambda,
            "_simulation.RDS"))
