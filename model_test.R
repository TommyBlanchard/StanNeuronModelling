library(rstan)

modelcode <- paste(readLines('multi_normal_covariances.stan'), collapse = '\n')

N <- 20

library(MASS) # For mvrnorm
data <- list(N=N,M=2,
             dim=2,
             alpha=c(1,1), # dirichlet prior
             x=mvrnorm(n=N, mu=c(0,0), Sigma=matrix(c(1,1.5, 1.5,3), nrow=2) )
             )


plot(data$x) # look at our data

myfit <- stan(model_code=modelcode, data=data, iter=1000, chains=1)


print(myfit)

plot(myfit)




