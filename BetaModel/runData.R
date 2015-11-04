library(rstan)
library(MASS) # For mvrnorm
source("../Shared.R")

filePath = '' #Set this to the path to the directory with your data

# make stan parallel 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

set.seed(100)
setup_data <- function(file_list, M){
d <- NULL
obs <- NULL
sigma <- array(0,dim=c(length(file_list),2,2))
for (file in 1:length(file_list)){
  celldata <-read.csv(file_list[file], header=FALSE);
  colnames(celldata) <- c('FR', 'X1', 'X2')
  fits <- lm(FR ~ X1 + X2, data = celldata);
  obs <- rbind(obs,fits$coefficients[2:3]);
  sigma[file,,] <- vcov(fits)[2:3,2:3]; #covariance matrix for the coefficients
}

alpha = 1; #Flat prior

# Construct the data to send
data <- list(N_CELLS=length(file_list), 
             M=M, 
             dim=2,
             alpha_cov_mix=array(rep(alpha,M),dim=M), # dirichlet prior. Need to declare as array to work with M=1
             alpha_noise_mix=c(1,1), # Set a uniform prior for signal vs noise
             x = obs,
             sigma = sigma
)      
return(data)
}

modelcode <- paste(readLines('beta_model.stan'), collapse = '\n')

file_list <- list.files(path=filePath, pattern="*.csv", full.names=TRUE)
runData = setup_data(file_list,1)
myfit <- stan(model_code=modelcode, data=runData, iter=1000, chains=4) #1000 iterations is probably the lower bound of what you should do - do higher if rhat is not below 1.1!
