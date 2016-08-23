library(rstan)
library(MASS) # For mvrnorm
source("../SharedNew.R")

set.seed(100)

make_data <- function(betas,N_CELLS,DATA_PER_CELL, noise){
  d <- NULL
  obs <- NULL
  sigma <- array(0,dim=c(N_CELLS,2,2))
  for(r in 1:nrow(betas)) {
    b <- betas[r,]
    
    m <- matrix(runif(DATA_PER_CELL*2), ncol=2)
    
    err <- rnorm(DATA_PER_CELL, 0.0, noise)
    d <- data.frame(cell=r, x1=m[,1], x2=m[,2], y= m%*%b+err)
    fits <- lm(scale(y) ~ scale(x1) + scale(x2), data = d);
    obs <- rbind(obs,fits$coefficients[2:3]);
    sigma[r,,] <- vcov(fits)[2:3,2:3]; #covariance matrix for the coefficients
  }
  M = 2
  print(M)
  # Construct the data to send
  data <- list(N_CELLS=N_CELLS, 
               M=M, 
               dim=2,
               alpha_cov_mix=array(rep(0.5,M),dim=M), # dirichlet prior. Need to declare as array to work with M=1
               X = obs,
               sigma = sigma
  )      
  return(data);
}

# Make up regression data

N_CELLS = 100

# 1 cov with positively correlated variables
betas_1poscov <- mvrnorm(n=N_CELLS, mu=c(0,0), Sigma=matrix(c(1,1.5, 1.5,3), nrow=2))
# 1 pos cov + noise (half and half)
betas_1poscov_noise <- rbind(mvrnorm(n=N_CELLS/2, mu=c(0,0), Sigma=matrix(c(1,1.5, 1.5,3), nrow=2)),mvrnorm(n=N_CELLS/2, mu=c(0,0), Sigma=matrix(c(0,0,0,0), nrow=2)))
# 1 circular cov
betas_1circcov <- mvrnorm(n=N_CELLS, mu=c(0,0), Sigma=matrix(c(1,0,0,1), nrow=2))

betas_1circcov <- mvrnorm(n=N_CELLS, mu=c(0,0), Sigma=matrix(c(3,0,0,3), nrow=2))
#x & y aligned
betas_xyaligned <- rbind(mvrnorm(n=N_CELLS/2, mu=c(0,0), Sigma=matrix(c(1,0, 0,0), nrow=2)),mvrnorm(n=N_CELLS/2, mu=c(0,0), Sigma=matrix(c(0,0, 0,1), nrow=2)))     
#x & y aligned + noise
betas_xyaligned_noise <- rbind(mvrnorm(n=N_CELLS/2, mu=c(0,0), Sigma=matrix(c(0,0, 0,0), nrow=2)),mvrnorm(n=N_CELLS/4, mu=c(0,0), Sigma=matrix(c(1,0, 0,0), nrow=2)),mvrnorm(n=N_CELLS/4, mu=c(0,0), Sigma=matrix(c(0,0, 0,1), nrow=2)))     

getMedianFit <- function(betas){
  # make stan parallel 
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores()) 
  
  modelcode <- paste(readLines('model.stan'), collapse = '\n')
  
  N_CELLS <- 100
  DATA_PER_CELL <- 100
  M <- 2;
  ITER <- 1000; # How many to run for
  CHAINS <- 4;
  
  data= make_data(betas,N_CELLS,DATA_PER_CELL, noise);
  fit <- stan(model_code=modelcode, data=data, iter=ITER, chains=CHAINS);
  
  
}


# data= make_data(mvrnorm(n=N_CELLS, mu=c(0,0), Sigma=matrix(c(0.1,0.1,0.1,.1), nrow=2)), N_CELLS, DATA_PER_CELL);
# fit <- stan(model_code=modelcode, data=data, iter=ITER, chains=CHAINS);


data= make_data(betas_1poscov_noise,N_CELLS,DATA_PER_CELL, 2);
fitposcov <- stan(model_code=modelcode, data=data, iter=ITER, chains=CHAINS);


data= make_data(betas_1circcov,N_CELLS,400, 1);
fit_corr <- stan(model_code=modelcode, data=data, iter=ITER, chains=CHAINS);

modelcode <- paste(readLines('model.stan'), collapse = '\n')
N_CELLS = 200
DATA_PER_CELL = 200
ITER = 1000
CHAINS = 4
 data= make_data(betas_xyaligned_noise,N_CELLS,DATA_PER_CELL, 2);
 fit2cov_2orthcov_noise <- stan(model_code=modelcode, data=data, iter=ITER, chains=CHAINS);
# 
# data= make_data(betas_2orthcov_noise,N_CELLS,DATA_PER_CELL);
 fit3cov_2orthcov_noise <- stan(model_code=modelcode, data=data, iter=ITER, chains=CHAINS);

 vis_dists(fitposcov, data, 'Variable 1', 'Variable 2', 'poscov')
 vis_posteriors(fitposcov, 'poscov') 
 