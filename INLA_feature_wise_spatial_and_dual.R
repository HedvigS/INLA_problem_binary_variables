
###

cat("#### Spatial only Model ####\n")

#make empty df to bind to
df_spatial_only <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(df_spatial_only) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "model", "waic", "marginals.hyperpar.spatial_id_iid_model", "marginals.hyperpar.spatial_id_generic") 
df_spatial_only$`2.5%` <- as.numeric(df_spatial_only$`2.5%`)
df_spatial_only$`50%` <- as.numeric(df_spatial_only$`50%`)
df_spatial_only$`97.5%` <- as.numeric(df_spatial_only$`97.5%`)
df_spatial_only$Feature_ID <- as.character(df_spatial_only$Feature_ID)
df_spatial_only$model <- as.character(df_spatial_only$model)
df_spatial_only$effect <- as.character(df_spatial_only$effect)
df_spatial_only$waic <- as.numeric(df_spatial_only$waic)
df_spatial_only$marginals.hyperpar.spatial_id_iid_model <- as.list(df_spatial_only$marginals.hyperpar.spatial_id_iid_model)
df_spatial_only$marginals.hyperpar.spatial_id_generic <- as.list(df_spatial_only$marginals.hyperpar.spatial_id_generic)


index <- 0

cat("Starting INLA spatial-only featurewise runs at", as.character(Sys.time()), ".\n")

cat("#### spatial only model ####\n")
for(feature in features){
  
  #feature <- features[1]
  
  cat(paste0("# Running the spatial-only model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n"))
  index <- index + 1 
  
  output <- try({eval(substitute(inla(formula = this_feature ~
                                        f((spatial_id_generic), 
                                          model = "generic0",
                                          Cmatrix = spatial_prec_mat,
                                          constr = TRUE, 
                                          hyper = pcprior) + 
                                        f(spatial_id_iid_model,
                                          model = "iid", 
                                          hyper = pcprior),
                                      control.compute = list(waic=TRUE, dic = FALSE, mlik = FALSE, config = TRUE),
                                      control.inla = list(tolerance = 1e-6, h = 0.001),
                                      control.predictor = list(compute=TRUE, link=1), #@Sam should we do this?
                                      control.family = list(control.link=list(model="logit")),   #@Sam should we do this?
                                      data = df,family = "binomial"),
                                 list(this_feature=as.name(feature))))   })
  
  if (class(output) != "try-error") {
    
    suppressWarnings(  saveRDS(output, file = paste0(OUTPUTDIR, "spatial_only/spatial_only_", feature, ".rdata")) )
    #Don't be alarmed by the suppress warnings. saveRDS() is being kind and reminding us that the package stats may not be available when loading. However, this is not a necessary warning for us so we've wrapped saveRDS in suppressWarnings.
    spatial_effect_generic = inla.tmarginal(function(x) 1/sqrt(x),
                                            output$marginals.hyperpar$`Precision for spatial_id_generic`,
                                            method = "linear") %>%
      inla.qmarginal(c(0.025, 0.5, 0.975), .)
    
    spatial_effect_iid_model = inla.tmarginal(function(x) 1/sqrt(x),
                                              output$marginals.hyperpar$`Precision for spatial_id_iid_model`,
                                              method = "linear") %>%
      inla.qmarginal(c(0.025, 0.5, 0.975), .)
    
    df_spatial_only_generic  <- spatial_effect_generic %>% 
      as.data.frame() %>% 
      t() %>% 
      as.data.frame() %>% 
      rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
      mutate(Feature_ID = feature) %>% 
      mutate(effect = "spatial_only_generic") %>% 
      mutate(model = "spatial_only") %>% 
      mutate(waic = output$waic$waic)  %>% 
      mutate(marginals.hyperpar.spatial_id_generic = output$marginals.hyperpar[1])
    
    df_spatial_only_iid_model <-   spatial_effect_iid_model %>% 
      as.data.frame() %>% 
      t() %>% 
      as.data.frame() %>% 
      rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
      mutate(Feature_ID = feature) %>% 
      mutate(effect = "spatial_only_iid_model") %>% 
      mutate(model = "spatial_only") %>% 
      mutate(waic = output$waic$waic)  %>% 
      mutate(marginals.hyperpar.spatial_id_iid_model = output$marginals.hyperpar[2])
    
    df_spatial_only <- df_spatial_only  %>% 
      full_join(df_spatial_only_iid_model, 
                by = c(join_columns, "marginals.hyperpar.spatial_id_iid_model")) %>% 
      full_join(df_spatial_only_generic, 
                by = c(join_columns, "marginals.hyperpar.spatial_id_generic"))
  }
  rm(output)  
}

df_spatial_only %>% write_tsv(file = file.path(OUTPUTDIR, "df_spatial_only.tsv"))
df_spatial_only %>% saveRDS(file = file.path(OUTPUTDIR, "df_spatial_only.Rdata"))

cat("All done with the spatial only model, 100% done!")

###
cat("#### AUTOTYP_area only model ####\n")

#make empty df to bind to
df_autotyp_area_only <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(df_autotyp_area_only) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "model", "waic", "marginals.hyperpar.AUTOTYP_area_id_iid_model") 
df_autotyp_area_only$`2.5%` <- as.numeric(df_autotyp_area_only$`2.5%`)
df_autotyp_area_only$`50%` <- as.numeric(df_autotyp_area_only$`50%`)
df_autotyp_area_only$`97.5%` <- as.numeric(df_autotyp_area_only$`97.5%`)
df_autotyp_area_only$Feature_ID <- as.character(df_autotyp_area_only$Feature_ID)
df_autotyp_area_only$model <- as.character(df_autotyp_area_only$model)
df_autotyp_area_only$effect <- as.character(df_autotyp_area_only$effect)
df_autotyp_area_only$waic <- as.numeric(df_autotyp_area_only$waic)
df_autotyp_area_only$marginals.hyperpar.AUTOTYP_area_id_iid_model <- as.list(df_autotyp_area_only$marginals.hyperpar.AUTOTYP_area_id_iid_model)


index <- 0

cat("Starting INLA AUTOTYP-area only runs featurewise runs at", as.character(Sys.time()), ".\n")

for(feature in features){
  
  #feature <- features[1]
  
  cat(paste0("# Running the autotyp_area-only model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n"))
  index <- index + 1 
  
  output <- try({eval(substitute(inla(formula = this_feature ~ 
                                        f(AUTOTYP_area_id_iid_model, 
                                          hyper = pcprior,
                                          model = "iid"),
                                      control.compute = list(waic=TRUE, dic = FALSE, mlik = FALSE, config = TRUE),
                                      control.predictor = list(compute=TRUE, link=1), #@Sam should we do this?
                                      control.family = list(control.link=list(model="logit")),   #@Sam should we do this?
                                      control.inla = list(tolerance = 1e-6, h = 0.001),
                                      data = df,family = "binomial"),
                                 list(this_feature=as.name(feature)))) })
  
  if (class(output) != "try-error") {
    
    suppressWarnings(    saveRDS(output, file = paste0(OUTPUTDIR, "autotyp_area_only/autotyp_area_only_", feature, ".rdata")))
    #Don't be alarmed by the suppress warnings. saveRDS() is being kind and reminding us that the package stats may not be available when loading. However, this is not a necessary warning for us so we've wrapped saveRDS in suppressWarnings
    
    autotyp_area_effect = inla.tmarginal(function(x) 1/sqrt(x),
                                         output$marginals.hyperpar$`Precision for AUTOTYP_area_id_iid_model`, method = "linear") %>%
      inla.qmarginal(c(0.025, 0.5, 0.975), .)
    
    df <- autotyp_area_effect %>% 
      as.data.frame() %>% 
      t() %>% 
      as.data.frame() %>% 
      rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
      mutate(Feature_ID = feature) %>% 
      mutate(effect = "autotyp_area_iid_model") %>% 
      mutate(model = "autotyp_area_only") %>% 
      mutate(waic = output$waic$waic)  %>% 
      mutate(marginals.hyperpar.AUTOTYP_area_id_iid_model = output$marginals.hyperpar[1])
    
    df_autotyp_area_only <- df_autotyp_area_only  %>% 
      full_join(df, by = c(join_columns, "marginals.hyperpar.AUTOTYP_area_id_iid_model"))
  }
  rm(output)  
}
cat("All done with the autotyp_area only model, 100% done!")

df_autotyp_area_only %>% write_tsv(file = file.path(OUTPUTDIR, "df_autotyp_area_only.tsv"))
df_autotyp_area_only %>% saveRDS(file = file.path(OUTPUTDIR, "df_autotyp_area_only.Rdata"))


cat("#### Spatial & Phylo Model ####\n")

index <- 0

#something was going awry with calculating the marginal effects of a specific features, GB051, so the for loop above has been split in twain: one which saves the entire output of inla() as an rdata object in a directory and one that calculates the effects and renders the same kind of df as above. This way running the for loop with inla() can still happen and we can debug the particulars after. 

cat("Starting INLA dual process (spatial and phylo) runs featurewise runs at", as.character(Sys.time()), ".\n")

for(feature in features){
  
  #feature <- features[2]
  
  cat(paste0("# Running the spatial-phylo (double-process) model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n"))
  index <- index + 1 
  
  output <- try({eval(substitute(inla(formula = this_feature ~
                                        f((phy_id_generic), 
                                          model = "generic0",
                                          Cmatrix = phylo_prec_mat,
                                          constr = TRUE, 
                                          hyper = pcprior) + 
                                        f(phy_id_iid_model,
                                          model = "iid", 
                                          hyper = pcprior,
                                          constr = TRUE) +
                                        f((spatial_id_generic), 
                                          model = "generic0",
                                          Cmatrix = spatial_prec_mat,
                                          constr = TRUE, 
                                          hyper = pcprior) + 
                                        f(spatial_id_iid_model,
                                          model = "iid", 
                                          hyper = pcprior,
                                          constr = TRUE),
                                      control.compute = list(waic=TRUE, dic = FALSE, mlik = FALSE, config = TRUE),
                                      control.inla = list(tolerance = 1e-6, h = 0.001),
                                      control.predictor = list(compute=TRUE, link=1), #@Sam should we do this?
                                      control.family = list(control.link=list(model="logit")),  #@Sam should we do this?
                                      data = df, family = "binomial"),
                                 list(this_feature=as.name(feature)))) })
  
  
  if (class(output) != "try-error") {
    suppressWarnings(    saveRDS(output, file = paste0(OUTPUTDIR, "dual_process_rdata/spatial_phylo_", feature, ".rdata")))
    #Don't be alarmed by the suppress warnings. saveRDS() is being kind and reminding us that the package stats may not be available when loading. However, this is not a necessary warning for us so we've wrapped saveRDS in suppressWarnings
    
  }
  
  rm(output)  
}

spatial_phylo_rdata_fns <- list.files(file.path(OUTPUTDIR, "/dual_process_rdata/"), full.names = T, pattern = ".*rdata")

#second for loop for the dual process model because of previously discussing debugging workflow the for loop is split in twain.

#make empty df to bind to
df_dual_spatial_phylo <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(df_dual_spatial_phylo) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "waic", "model", "marginals.hyperpar.phy_id_iid_model", "marginals.hyperpar.phy_id_generic", "marginals.hyperpar.spatial_id_iid_model", "marginals.hyperpar.spatial_id_generic") 
df_dual_spatial_phylo$`2.5%` <- as.numeric(df_dual_spatial_phylo$`2.5%`)
df_dual_spatial_phylo$`50%` <- as.numeric(df_dual_spatial_phylo$`50%`)
df_dual_spatial_phylo$Feature_ID <- as.character(df_dual_spatial_phylo$Feature_ID)
df_dual_spatial_phylo$effect <- as.character(df_dual_spatial_phylo$effect)
df_dual_spatial_phylo$model <- as.character(df_dual_spatial_phylo$model)
df_dual_spatial_phylo$waic <- as.numeric(df_dual_spatial_phylo$waic)
df_dual_spatial_phylo$marginals.hyperpar.phy_id_iid_model <- as.list(df_dual_spatial_phylo$marginals.hyperpar.phy_id_iid_model)
df_dual_spatial_phylo$marginals.hyperpar.phy_id_generic <- as.list(df_dual_spatial_phylo$marginals.hyperpar.phy_id_generic)
df_dual_spatial_phylo$marginals.hyperpar.spatial_id_iid_model <- as.list(df_dual_spatial_phylo$marginals.hyperpar.spatial_id_iid_model)
df_dual_spatial_phylo$marginals.hyperpar.spatial_id_generic <- as.list(df_dual_spatial_phylo$marginals.hyperpar.spatial_id_generic)

for(fn in spatial_phylo_rdata_fns){
  
  #  fn <- spatial_phylo_rdata_fns[1]
  
  output <- readRDS(fn)
  
  feature <- fn %>% str_extract("GB[0-9]*[a|b]?")
  
  cat(paste0("I'm processing the inla output for feature ", feature, ".\n" ))
  
  phylo_effect_generic = inla.tmarginal(function(x) 1/sqrt(x),
                                        output$marginals.hyperpar$`Precision for phy_id_generic`,
                                        method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  phylo_effect_iid_model = inla.tmarginal(function(x) 1/sqrt(x),
                                          output$marginals.hyperpar$`Precision for phy_id_iid_model`,
                                          method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  
  spatial_effect_generic = inla.tmarginal(function(x) 1/sqrt(x),
                                          output$marginals.hyperpar$`Precision for spatial_id_generic`,
                                          method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  spatial_effect_iid_model = inla.tmarginal(function(x) 1/sqrt(x),
                                            output$marginals.hyperpar$`Precision for spatial_id_iid_model`,
                                            method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  df_phylo_generic  <- phylo_effect_generic %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "phylo_generic_in_dual") %>% 
    mutate(model = "dual") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.phy_id_generic = output$marginals.hyperpar[1])
  
  df_phylo_iid_model <-   phylo_effect_iid_model %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "phylo_iid_model_in_dual") %>% 
    mutate(model = "dual") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.phy_id_iid_model = output$marginals.hyperpar[2])
  
  df_spatial_generic  <- spatial_effect_generic %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "spatial_generic_in_dual") %>% 
    mutate(model = "dual") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.spatial_id_generic = output$marginals.hyperpar[3])
  
  df_spatial_iid_model <-   spatial_effect_iid_model %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "spatial_iid_model_in_dual") %>% 
    mutate(model = "dual") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.spatial_id_iid_model = output$marginals.hyperpar[4])
  
  df_dual_spatial_phylo <- df_dual_spatial_phylo %>%
    full_join(df_phylo_generic, by = c(join_columns,
                                       "marginals.hyperpar.phy_id_generic")) %>% 
    full_join(df_phylo_iid_model, by = c(join_columns,
                                         "marginals.hyperpar.phy_id_iid_model")) %>% 
    full_join(df_spatial_generic, by = c(join_columns, 
                                         "marginals.hyperpar.spatial_id_generic")) %>% 
    full_join(df_spatial_iid_model, by = c(join_columns,
                                           "marginals.hyperpar.spatial_id_iid_model"))
  rm(output)  
}

df_dual_spatial_phylo %>% write_tsv(file = file.path(OUTPUTDIR, "df_spatial_phylo.tsv"))
df_dual_spatial_phylo %>% saveRDS(file = file.path(OUTPUTDIR, "df_spatial_phylo.Rdata"))

#TRIAL PROCESS

cat("#### Spatial & Phylo Model + AUTOTYP area ####\n")

index <- 0

cat("Starting INLA trial process (AUTOTYP-area, spatial and phylo) runs featurewise runs at", as.character(Sys.time()), ".\n")

for(feature in features){
  
  #feature <- features[2]
  
  cat(paste0("# Running the spatial-phylo-area (trial-process) model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n"))
  index <- index + 1 
  
  output <- try({eval(substitute(inla(formula = this_feature ~
                                        f((phy_id_generic), 
                                          model = "generic0",
                                          Cmatrix = phylo_prec_mat,
                                          constr = TRUE, 
                                          hyper = pcprior) + 
                                        f(phy_id_iid_model,
                                          model = "iid", 
                                          hyper = pcprior) +
                                        f((spatial_id_generic), 
                                          model = "generic0",
                                          Cmatrix = spatial_prec_mat,
                                          constr = TRUE, 
                                          hyper = pcprior) + 
                                        f(spatial_id_iid_model,
                                          model = "iid", 
                                          hyper = pcprior) +  
                                        f(AUTOTYP_area_id_iid_model,
                                          model = "iid",
                                          hyper = pcprior),
                                      control.compute = list(waic=TRUE, dic = FALSE, mlik = FALSE, config = TRUE),
                                      control.predictor = list(compute=TRUE, link=1), #@Sam should we do this?
                                      control.family = list(control.link=list(model="logit")),   #@Sam should we do this?
                                      control.inla = list(tolerance = 1e-6, h = 0.001),
                                      data = df, family = "binomial"),
                                 list(this_feature=as.name(feature)))) })
  
  if (class(output) != "try-error") {
    
    suppressWarnings(saveRDS(output, file = paste0(OUTPUTDIR, "/trial_process_rdata/spatial_phylo_area_", feature, ".rdata")))
    #Don't be alarmed by the suppress warnings. saveRDS() is being kind and reminding us that the package stats may not be available when loading. However, this is not a necessary warning for us so we've wrapped saveRDS in suppressWarnings
  }
  rm(output)  
}

spatial_phylo_area_rdata_fns <- list.files(file.path(OUTPUTDIR, "trial_process_rdata/"), full.names = T, pattern = ".*rdata")

#second for loop for the dual process model because of previously discussing debugging workflow the for loop is split in twain.

#make empty df to bind to
df_trial <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(df_trial) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "waic", "model",  "marginals.hyperpar.phy_id_iid_model", "marginals.hyperpar.phy_id_generic", "marginals.hyperpar.spatial_id_iid_model", "marginals.hyperpar.spatial_id_generic", "marginals.hyperpar.AUTOTYP_area_id_iid_model") 
df_trial $`2.5%` <- as.numeric(df_trial $`2.5%`)
df_trial $`50%` <- as.numeric(df_trial $`50%`)
df_trial $`97.5%` <- as.numeric(df_trial $`97.5%`)
df_trial $Feature_ID <- as.character(df_trial $Feature_ID)
df_trial $effect <- as.character(df_trial $effect)
df_trial $waic <- as.numeric(df_trial$waic)
df_trial $model <- as.character(df_trial$model)
df_trial$marginals.hyperpar.phy_id_iid_model <- as.list(df_trial$marginals.hyperpar.phy_id_iid_model)
df_trial$marginals.hyperpar.phy_id_generic <- as.list(df_trial$marginals.hyperpar.phy_id_generic)
df_trial$marginals.hyperpar.spatial_id_iid_model <- as.list(df_trial$marginals.hyperpar.spatial_id_iid_model)
df_trial$marginals.hyperpar.spatial_id_generic <- as.list(df_trial$marginals.hyperpar.spatial_id_generic)
df_trial$marginals.hyperpar.AUTOTYP_area_id_iid_model <- as.list(df_trial$marginals.hyperpar.AUTOTYP_area_id_iid_model)

for(fn in spatial_phylo_area_rdata_fns){
  
  #  fn <- spatial_phylo_area_rdata_fns[1]
  # fn <- "output/spatiophylogenetic_modelling/results_small/dual_process_rdata/spatial_phylo_area_GB051.rdata"
  
  output <- readRDS(fn)
  
  feature <- fn %>% str_extract("GB[0-9]*[a|b]?")
  
  cat(paste0("I'm processing the inla output for feature ", feature, "in the trial model.\n" ))
  
  
  phylo_effect_generic = inla.tmarginal(function(x) 1/sqrt(x),
                                        output$marginals.hyperpar$`Precision for phy_id_generic`,
                                        method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  phylo_effect_iid_model = inla.tmarginal(function(x) 1/sqrt(x),
                                          output$marginals.hyperpar$`Precision for phy_id_iid_model`,
                                          method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  
  spatial_effect_generic = inla.tmarginal(function(x) 1/sqrt(x),
                                          output$marginals.hyperpar$`Precision for spatial_id_generic`,
                                          method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  spatial_effect_iid_model = inla.tmarginal(function(x) 1/sqrt(x),
                                            output$marginals.hyperpar$`Precision for spatial_id_iid_model`,
                                            method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  autotyp_area_effect_iid_model = inla.tmarginal(function(x) 1/sqrt(x),
                                                 output$marginals.hyperpar$`Precision for AUTOTYP_area_id_iid_model`,
                                                 method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  
  df_phylo_generic  <- phylo_effect_generic %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "phylo_generic_in_trial") %>% 
    mutate(model = "trial") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.phy_id_generic = output$marginals.hyperpar[1])
  
  df_phylo_iid_model <-   phylo_effect_iid_model %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "phylo_iid_model_in_trial") %>% 
    mutate(model = "trial") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.phy_id_iid_model = output$marginals.hyperpar[2])
  
  df_spatial_generic  <- spatial_effect_generic %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "spatia_generic_in_trial") %>% 
    mutate(model = "trial") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.spatial_id_generic = output$marginals.hyperpar[3])
  
  df_spatial_iid_model <-   spatial_effect_iid_model %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "spatial_iid_model_in_trial") %>% 
    mutate(model = "trial") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.spatial_id_iid_model = output$marginals.hyperpar[4])
  
  df_autotyp_iid  <- autotyp_area_effect_iid_model %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "autotyp_iid_in_trial") %>% 
    mutate(model = "trial") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.AUTOTYP_area_id_iid_model = output$marginals.hyperpar[5])
  
  df_trial <- df_trial  %>% 
    full_join(df_phylo_generic, 
              by = c(join_columns ,"marginals.hyperpar.phy_id_generic")) %>% 
    full_join(df_phylo_iid_model, 
              by = c(join_columns,"marginals.hyperpar.phy_id_iid_model")) %>% 
    full_join(df_spatial_generic, 
              by = c(join_columns, "marginals.hyperpar.spatial_id_generic")) %>% 
    full_join(df_spatial_iid_model, 
              by = c(join_columns,"marginals.hyperpar.spatial_id_iid_model")) %>% 
    full_join(df_autotyp_iid, 
              by = c(join_columns,"marginals.hyperpar.AUTOTYP_area_id_iid_model"))
  
  
  rm(output)  
}

df_trial %>% write_tsv(file = file.path(OUTPUTDIR, "df_trial.tsv"))
df_trial %>% saveRDS(file = file.path(OUTPUTDIR, "df_trial.Rdata"))

cat("All done with all INLA runs at ", as.character(Sys.time())[1], ".\n")

