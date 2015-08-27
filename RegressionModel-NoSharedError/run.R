library(rstan)
library(MASS) # For mvrnorm

source("../Shared.R")

# make stan parallel 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

options(width=Sys.getenv("COLUMNS"))  # fix the number of columns

modelcode <- paste(readLines('model.stan'), collapse = '\n')

N_CELLS <- 100
DATA_PER_CELL <- 100

s <- rep.int(DATA_PER_CELL,N_CELLS) # Number of responses in each cell

set.seed(100)

# Make up regression data

# 1 cov with positively correlated variables
betas <- mvrnorm(n=N_CELLS, mu=c(0,0), Sigma=matrix(c(1,1.5, 1.5,3), nrow=2))

# 2 covs, orthogonal to each other
betas <- rbind(mvrnorm(n=N_CELLS/2, mu=c(0,0), Sigma=matrix(c(1,1.5, 1.5,3), nrow=2)),mvrnorm(n=N_CELLS/2, mu=c(0,0), Sigma=matrix(c(1,-1.5, -1.5,3), nrow=2)))

# 1 cov + noise (half and half)
betas <- rbind(mvrnorm(n=N_CELLS/2, mu=c(0,0), Sigma=matrix(c(1,1.5, 1.5,3), nrow=2)),mvrnorm(n=N_CELLS/2, mu=c(0,0), Sigma=matrix(c(0,0,0,0), nrow=2)))

# 1 circular cov
betas <- mvrnorm(n=N_CELLS, mu=c(0,0), Sigma=matrix(c(3,0,0,3), nrow=2))

# 2 orthogonal covs + noise (half noise)
betas <- rbind(mvrnorm(n=N_CELLS/4, mu=c(0,0), Sigma=matrix(c(1,1.5, 1.5,3), nrow=2)),mvrnorm(n=N_CELLS/4, mu=c(0,0), Sigma=matrix(c(1,-1.5, -1.5,3), nrow=2)),mvrnorm(n=N_CELLS/2, mu=c(0,0), Sigma=matrix(c(0,0,0,0), nrow=2)))


# Show our betas
plot(betas)
abline(h=0,v=0)

d <- NULL
for(r in 1:nrow(betas)) {
    b <- betas[r,]
    
    m <- matrix(runif(DATA_PER_CELL*2), ncol=2)
    
    err <- rnorm(DATA_PER_CELL, 0.0, 0.10)
    d <- rbind(d, data.frame(cell=r, x1=m[,1], x2=m[,2], y= m%*%b+err)) 
}

# If we wanted to recover the coefficients
# l <- lm( y ~ (x1 + x2) * as.factor(cell) -x1-x2, data=d)

# Construct the data to send
data <- list(N_CELLS=N_CELLS, 
             s=s,
             M=3, 
             N_RESPONSES=nrow(d),
             dim=2,
             alpha_cov_mix=c(1,1,1)*.1, # dirichlet prior
             alpha_noise_mix=c(1,10), # Set a prior favoring signal over noise
             x1=d$x1, x2=d$x2, cell=d$cell, y=d$y
)              


# plot(data$beta) # look at our data

myfit <- stan(model_code=modelcode, data=data, iter=1000, chains=1) #, control=list(stepsize=0.001))
# myfit <- stan(model_code=modelcode, data=data, iter=10000, warmup=10, chains=1, control=list(adapt_engaged=FALSE), algorithm="HMC")
# myfit <- stan(model_code=modelcode, data=data, iter=10000, warmup=10, chains=1, algorithm="HMC")
# myfit <- stan(model_code=modelcode, data=data, iter=500000, warmup=5, chains=1) # Holy crap warmup takes forever

# print(myfit)
# 
 plot(myfit)

# traceplot(myfit, pars = 'beta')

# traceplot(myfit, pars = 'residual_variance')

print(myfit, pars='beta')

print(myfit, pars='covs')

print(myfit, pars='mixture_weights')

print(myfit,pars='noise_weight')

print(myfit, pars='scale')

print(myfit, pars='noise_error')
print(myfit, pars='signal_error')

plot_beta_dif(myfit, betas)
# plot_covs(myfit, 4)

# get_sampler_params(myfit)


summary(do.call(rbind, args = get_sampler_params(myfit, inc_warmup = TRUE)), digits = 2)