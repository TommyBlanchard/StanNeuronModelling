
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
    vector[2] alpha_noise_mix; //beta prior
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
   
   real<lower=0> signal_error;  // half cauchy
   real<lower=0> noise_error;  // half cauchy
 
   vector<lower=0>[2] scale; // half normal
}
    
model {    
    real logp_cov_beta[M]; // for saving log probs
  
    real logps[M]; // for saving log probs
    real logpsnoise[2]; // noise or not?

    vector[N_CELLS] logp_noise; // for saving log prob of responses given neuron is noise
    vector[N_CELLS] logp_beta;  // for saving log prob of responses given a neuron's betas

    // P(signal_error)
    signal_error ~ cauchy(0, 1);
    noise_error  ~ cauchy(0, 1);
    
    // P(scaling factors)
    scale ~ cauchy(0, 1);
    
    // P(mixture_weights)
    mixture_weights ~ dirichlet(alpha_cov_mix);
    
    // P(noise_weight)
    noise_weight ~ beta(alpha_noise_mix[1], alpha_noise_mix[2]);
    
    // P(cov)
    for (m in 1:M) {
        covs[m] ~ inv_wishart(dim, identity);
    }
    
    // P(beta_c | covs) -- marginalize over covs for computing the betas
    for (n in 1:N_CELLS) {
        
        beta0[n] ~ normal(0,1);
        
        // but then penalize by summing over covariances
        // Directly compute the marginalized log likelihood, since stan can't handle discrete variables
        // see: http://www.michaelchughes.com/blog/2012/09/review-of-stan-off-the-shelf-hamiltonian-mcmc/
        for (m in 1:M) {
           // prob of the betas
           logp_cov_beta[m] <- log(mixture_weights[m]) + multi_normal_log(beta[n], zeros, diag_matrix(scale)*covs[m]*diag_matrix(scale));
        }
        increment_log_prob(log_sum_exp(logp_cov_beta));
    }
    
    logp_noise <- rep_vector(0.0, N_CELLS);
    logp_beta  <- rep_vector(0.0, N_CELLS); 
    
    // P( y_n | beta_c) -- given those betas, how likely is the data?
    // well this comes from marginalizing over whether or not each cell is noise (meaning we must take the product of its responses)
    for(r in 1:N_RESPONSES) {
        logp_noise[cell[r]] <- logp_noise[cell[r]] + normal_log(y[r], 0, noise_error);
        logp_beta[cell[r]]  <- logp_beta[cell[r]]  + normal_log(y[r], beta0[cell[r]] + beta[cell[r],1]*x1[r] + beta[cell[r],2]*x2[r], signal_error); 
        // y[r] ~ normal( beta0[cell[r]] + beta[cell[r],1]*x1[r] + beta[cell[r],2]*x2[r], signal_error);
    }

    for (n in 1:N_CELLS) {   
        increment_log_prob(log_sum_exp(log(noise_weight) + logp_noise[n], log(1.0-noise_weight) + logp_beta[n]));
    }
    
    //print("covs=", covs); 
}
