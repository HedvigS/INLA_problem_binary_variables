## Run Script

## Signal strength parameters
lambda_strong = 0.8
lambda_medium = 0.5
lambda_weak = 0.2

CLI <- "no" #set to "Yes" if you want it
brms <- "no" #set to "Yes" if you want it

if(CLI == "Yes") {
  command1 = paste0("RScript single_lambda_test.R ", lambda_strong)
  command2 = paste0("RScript single_geography_test.R ", lambda_strong)
  command3 = paste0("RScript dual_test.R ", lambda_strong)
  
  system(command1)
  system(command2)
  system(command3)
  
  ## Medium Signal
  command4 = paste0("RScript single_lambda_test.R ", lambda_medium)
  command5 = paste0("RScript single_geography_test.R ", lambda_medium)
  command6 = paste0("RScript dual_test.R ", lambda_medium)
  
  system(command4)
  system(command5)
  system(command6)
  
  ## Weak Signal
  command7 = paste0("RScript single_lambda_test.R ", lambda_weak)
  command8 = paste0("RScript single_geography_test.R ", lambda_weak)
  command9 = paste0("RScript dual_test.R ", lambda_weak)
  
  system(command7)
  system(command8)
  system(command9)
  
} else {
  
  source("install_inla.R")
  
  sink(file = "output/log_file_simulation_test.txt", split = T)
  lambda <- lambda_strong
  source("simulation_test/single_lambda_test.R")
  source("simulation_test/single_geography_test.R")
  source("simulation_test/dual_test.R")
    
  lambda <- lambda_medium
  source("simulation_test/single_lambda_test.R")
  source("simulation_test/single_geography_test.R")
  source("simulation_test/dual_test.R")
  
  lambda <- lambda_weak
  source("simulation_test/single_lambda_test.R")
  source("simulation_test/single_geography_test.R")
  source("simulation_test/dual_test.R")
  
  sink(file = NULL)
  
}

