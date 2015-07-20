library(rstan)
library(MASS) # For mvrnorm
options(width=Sys.getenv("COLUMNS")) 
 
modelcode <- paste(readLines('model.stan'), collapse = '\n')

N_CELLS <- 50
DATA_PER_CELL <- 150

# Make up regression data
betas <- mvrnorm(n=N_CELLS, mu=c(0,0), Sigma=matrix(c(1,1.5, 1.5,3), nrow=2))

d <- NULL
for(r in 1:nrow(betas)) {
    b <- betas[r,]
    
    m <- matrix(runif(DATA_PER_CELL), ncol=2)
    
    d <- rbind(d, data.frame(cell=r, x1=m[,1], x2=m[,2], y=(m %*% b)))
}

# Construct the data to send
data <- list(N_CELLS=N_CELLS, M=2,
             N_RESPONSES=nrow(d),
             dim=2,
             alpha=c(1,1), # dirichlet prior
             x1=d$x1, x2=d$x2, cell=d$cell, y=d$y
             )            


# plot(data$beta) # look at our data

myfit <- stan(model_code=modelcode, data=data, iter=10000, chains=4)


# print(myfit)
# 
# plot(myfit)

# traceplot(myfit, pars = 'beta')

traceplot(myfit, pars = 'residual_variance')

print(myfit, pars='betas')
print(myfit, pars='covs')

