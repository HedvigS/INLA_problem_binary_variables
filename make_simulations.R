## Run Script

## Strong Signal
lambda_strong = 0.8

command1 = paste0("RScript single_lambda_test.R ", lambda_strong)
command2 = paste0("RScript single_geography_test.R ", lambda_strong)

system(command1)
system(command2)

## Medium Signal
lambda_medium = 0.5

command3 = paste0("RScript single_lambda_test.R ", lambda_medium)
command4 = paste0("RScript single_geography_test.R ", lambda_medium)

system(command4)
system(command5)

## Weak Signal
lambda_weak = 0.2

command7 = paste0("RScript single_lambda_test.R ", lambda_weak)
command8 = paste0("RScript single_geography_test.R ", lambda_weak)

system(command3)
system(command4)