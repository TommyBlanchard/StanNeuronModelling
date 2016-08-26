library(rstan)
library(MASS) # For mvrnorm
library(plyr)
library(coda)

dataDir = "../singleunitdata/";
load(file=paste(dataDir, 'regdata.RData', sep = ""))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

readyData <- function(dat){
  M = 2
  alpha = 0.5
  list(
    N_CELLS = length(dat$neuron),
    M = M,
    dim = 2,
    alpha_cov_mix = array(rep(alpha,M),dim = M), # dirichlet prior. Need to declare as array to work with M=1
    alpha_noise_mix = c(0.5,0.5),
    X = cbind(dat$beta1,dat$beta2),
    sigma = array(c(dat$sigma11,dat$sigma21,dat$sigma12,dat$sigma22), dim=c(length(dat$sigma11),2,2))
  )
}

fitStanModel <- function(dat){
  set.seed(100)
  ITER <- 5000; #1000 iterations is probably the lower bound of what you should do - do higher if rhat is not below 1.1!
  CHAINS <- 5;
  modelcode <- paste(readLines('../Model/model.stan'), collapse = '\n')
  
  alpha = 0.5; #alpha for dirichlet prior
  
  # Construct the data to send
  data <- readyData(dat)
  
  c(data, stan(model_code=modelcode, data=data, iter=ITER, chains=CHAINS))
}

fits = dlply(regdata,'setName',fitStanModel)

save(fits,file=paste(dataDir, 'fiveKfits.RData', sep = ""))

