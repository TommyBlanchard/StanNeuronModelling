// List of matrices in R: https://groups.google.com/forum/#!topic/stan-users/-6s9PvVsyug

data {
    int<lower=0> N; // Sample size
    int<lower=1> dim; // Number of dimensions
    int<lower=1> M; // number of mixture components 
    row_vector[dim] x[N]; // Value for each sample on each dimension
    
    vector[M] alpha; // dirichlet prior
}

transformed data {
    vector[dim] zeros; //For now, means are zero
    matrix[dim,dim] identity; //Identity for scaling matrix
    
    zeros <- rep_vector(0, dim);
    identity <- diag_matrix(rep_vector(1.0,dim)); 
}

parameters {
   cov_matrix[dim] covs[M];
   
   simplex[M] mixture_weights; 
}
    
model {    
    real logps[M]; // for saving log probs
    
    mixture_weights ~ dirichlet(alpha);
        
    for (m in 1:M)
        covs[m] ~ inv_wishart(dim, identity);
        
    // Directly compute the marginalized log likelihood, since stan can't handle discrete variables
    // see: http://www.michaelchughes.com/blog/2012/09/review-of-stan-off-the-shelf-hamiltonian-mcmc/
    for (n in 1:N) {
        for (m in 1:M) {
            logps[m] <- log(mixture_weights[m]) + multi_normal_log(x[n], zeros, covs[m]);
        }
        lp__ <- lp__ + log_sum_exp(logps);
    }
    
    //print("covs=", covs); 
}

