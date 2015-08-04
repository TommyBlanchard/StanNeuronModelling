
data {
    int<lower=0> N_CELLS; // Sample size
    int<lower=0> N_RESPONSES; // number of total responses (N_CELLS * responses per)
    int<lower=1> dim; // Number of dimensions
    int<lower=1> M; // number of mixture components. The first is constrained to be essentially zero mean
    
    int  cell[N_RESPONSES]; // which cell is each response (grouping factor; MUST START AT 1 AND INCREMENT BY 1)
 
    real   x1[N_RESPONSES]; // predictors x1, x2
    real   x2[N_RESPONSES]; // predictors x1, x2
    real    y[N_RESPONSES]; // the response
    
    vector[M] alpha_cov_mix; // dirichlet priors
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
   real<lower=0, upper=1> noise_weight;
    
   real           beta0[N_CELLS]; // intercepts, not modeled with covariance. HMM Should they be?
   row_vector[dim] beta[N_CELLS];
   real<lower=0>  noise[N_CELLS];
   
   vector<lower=0>[2] scale; // half normal
}
    
model {    
    real logp_cov_beta[M]; // for saving log probs
  
    real logps[M]; // for saving log probs
    real logpsnoise[2]; // noise or not?

    // P(scaling factors)
    scale ~ cauchy(0, 10);
    
    // P(mixture_weights)
    mixture_weights ~ dirichlet(alpha_cov_mix);
    
    // P(cov)
    for (m in 1:M) {
        covs[m] ~ inv_wishart(dim, identity);
    }
    
    // P(beta_c | covs) -- marginalize over covs for computing the betas
    for (n in 1:N_CELLS) {
        
        beta0[n] ~ normal(0,10);
        noise[n] ~ cauchy(0,10);
        
        // but then penalize by summing over covariances
        // Directly compute the marginalized log likelihood, since stan can't handle discrete variables
        // see: http://www.michaelchughes.com/blog/2012/09/review-of-stan-off-the-shelf-hamiltonian-mcmc/
        for (m in 1:M) {
           // prob of the betas
           logp_cov_beta[m] <- log(mixture_weights[m]) + multi_normal_log(beta[n], zeros, diag_matrix(scale)*covs[m]*diag_matrix(scale));
        }
        increment_log_prob(log_sum_exp(logp_cov_beta));
    }
    
    // P( y_n | beta_c) -- given those betas, how likely is the data?
    // well this comes from marginalizing over whether or not each cell is noise (meaning we must take the product of its responses)
    for(r in 1:N_RESPONSES) {
        y[r] ~ normal( beta0[cell[r]] + beta[cell[r],1]*x1[r] + beta[cell[r],2]*x2[r], noise[cell[r]]);
    }

    //print("covs=", covs); 
}
