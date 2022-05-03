## Run Script

## Strong Signal
lambda_strong = 0.8
lambda_medium = 0.5
lambda_weak = 0.2

CLI <- "no"

if(CLI == "Yes") {
command1 = paste0("RScript simulation_test/single_lambda_test.R ", lambda_strong)
command2 = paste0("RScript simulation_test/single_geography_test.R ", lambda_strong)

system(command1)
system(command2)

## Medium Signal
command3 = paste0("RScript simulation_test/single_lambda_test.R ", lambda_medium)
command4 = paste0("RScript simulation_test/single_geography_test.R ", lambda_medium)

system(command4)
system(command5)

## Weak Signal
command7 = paste0("RScript simulation_test/single_lambda_test.R ", lambda_weak)
command8 = paste0("RScript simulation_test/single_geography_test.R ", lambda_weak)

system(command3)
system(command4)

} else{
  
source("install_inla.R")
  
sink(file = "output/log_file_simulation_test.txt", split = T)
  lambda <- lambda_strong
  source("simulation_test/single_lambda_test.R")
  source("simulation_test/single_geography_test.R")

lambda <- lambda_medium
  source("simulation_test/single_lambda_test.R")
  source("simulation_test/single_geography_test.R")

lambda <- lambda_weak
source("simulation_test/single_lambda_test.R")
source("simulation_test/single_geography_test.R")
  
sink(file = NULL)
    
  }

