library(rstan)
library(MASS) # For mvrnorm
library(plyr)
library(coda)

args <- commandArgs(trailingOnly = TRUE)

startNeurons = args[1]
betaType = args[2]

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


fitWeights <- function(betas, neurons, compileFit){
  
  modelcode <- paste(readLines('model.stan'), collapse = '\n')
  
  M <- 2;
  ITER <- 1000; # How many to run for
  CHAINS <- 4;
  noise = .3;
  DATA_PER_CELL = 500;
  data= make_data(betas,neurons,DATA_PER_CELL, noise);
  fit <- stan(model_code=modelcode, data=data, fit=compileFit, iter=ITER, chains=CHAINS);
  
  pars = extract(fit, permuted='true')
  
  free_weight = pars$mix_weights[,2]
  free_weight_med = median(free_weight)
  free_weight_high= HPDinterval(as.mcmc(free_weight), prob=0.95)[2];
  free_weight_low = HPDinterval(as.mcmc(free_weight), prob=0.95)[1];
  
  noise_weight = array(pars$noise_weight)
  noise_weight_med = median(noise_weight)
  noise_weight_high= HPDinterval(as.mcmc(noise_weight), prob=0.95)[2];
  noise_weight_low = HPDinterval(as.mcmc(noise_weight), prob=0.95)[1];
  cbind(free_weight_med, free_weight_low, free_weight_high, noise_weight_med, noise_weight_low, noise_weight_high, neurons)
}

set.seed(100)

if(file.exists(paste('weights', startNeurons, '.', betaType, '.RData', sep=''))){
	load(paste('weights', startNeurons, '.', betaType, '.RData', sep=''))
}else{
	i = 1
	d = NULL
}

if (startNeurons == 1){
	i = 10;
}

if (i < 20){
betas <- mvrnorm(n=10, mu=c(0,0), Sigma=matrix(c(1,0,0,1),nrow=2))
data <- make_data(betas,10,10,0) 
  modelcode <- paste(readLines('model.stan'), collapse = '\n')
fit <- stan(model_code=modelcode, data=data, iter=1, chains=1);
for (i in i:20){
  N_CELLS <- (as.integer(startNeurons)-1)*20 + i
  print(N_CELLS)
  print(i)
  print(d)
save(d, i, file=paste('weights', startNeurons, '.', betaType, '.RData', sep=''))

  if (betaType == 0){
    # 1 cov with positively correlated variables
    betas <- mvrnorm(n=N_CELLS, mu=c(0,0), Sigma=matrix(c(1,0.5, 0.5,1), nrow=2))
  }else if (betaType == 1){
    # 1 pos cov + noise (half and half)
    betas <- rbind(mvrnorm(n=ceiling(N_CELLS/2), mu=c(0,0), Sigma=matrix(c(.1,0, 0,.1), nrow=2)),mvrnorm(n=floor(N_CELLS/2), mu=c(0,0), Sigma=matrix(c(0,0,0,0), nrow=2)))
  }else if (betaType == 2){
    #x & y aligned + noise
    betas <- rbind(mvrnorm(n=ceiling(N_CELLS/2), mu=c(0,0), Sigma=matrix(c(0,0, 0,0), nrow=2)),mvrnorm(n=ceiling(floor(N_CELLS/2)/2), mu=c(0,0), Sigma=matrix(c(.1,0, 0,0), nrow=2)),mvrnorm(n=floor(floor(N_CELLS/2)/2), mu=c(0,0), Sigma=matrix(c(0,0, 0,.1), nrow=2)))     
  }else if (betaType == 3){
    #x & y aligned + circ cov + noise
    betas <- rbind(mvrnorm(n=floor(N_CELLS/4), mu=c(0,0), Sigma=matrix(c(.1,0, 0,.1), nrow=2)),mvrnorm(n=floor(N_CELLS/4), mu=c(0,0), Sigma=matrix(c(.1,0, 0,0), nrow=2)),mvrnorm(n=floor(N_CELLS/4), mu=c(0,0), Sigma=matrix(c(0,0, 0,.1), nrow=2)), mvrnorm(n=floor(N_CELLS/4) + N_CELLS%%4, mu=c(0,0), Sigma=matrix(c(0,0, 0,0), nrow=2)))     
  }
  d = rbind(d, fitWeights(betas, N_CELLS, fit))
  }
}

