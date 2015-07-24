// List of matrices in R: https://groups.google.com/forum/#!topic/stan-users/-6s9PvVsyug

data {
    int<lower=0> N_CELLS; // Sample size
    int<lower=0> N_RESPONSES; // number of total responses (N_CELLS * responses per)
    int<lower=1> dim; // Number of dimensions
    int<lower=1> M; // number of mixture components. The first is constrained to be essentially zero mean
    
<<<<<<< HEAD
    int  cell[N_RESPONSES]; // which cell is each response (grouping factor; MUST START AT 1 AND INCREMENT BY 1)
=======
    int  cell[N_RESPONSES]; // which cell is each response (grouping factor)
>>>>>>> c5eb4f20bc0b8ebf231fad876f98e7645cba8635
    real   x1[N_RESPONSES]; // predictors x1, x2
    real   x2[N_RESPONSES]; // predictors x1, x2
    real    y[N_RESPONSES]; // the response
    
<<<<<<< HEAD
    vector[M] alpha_cov_mix; // dirichlet priors
    vector[2] alpha_noise_mix; //beta prior
=======
    vector[M] alpha; // dirichlet prior
>>>>>>> c5eb4f20bc0b8ebf231fad876f98e7645cba8635
}

transformed data {
    vector[dim] zeros; //For now, means are zero
    matrix[dim,dim] identity; //Identity for scaling matrix
   
    zeros <- rep_vector(0, dim);
    identity <- diag_matrix(rep_vector(100.0,dim));    
}

parameters {
   cov_matrix[dim] covs[M];
   
   simplex[M] mixture_weights; 
<<<<<<< HEAD
   real<lower=0, upper=1> noise_weight;
=======
>>>>>>> c5eb4f20bc0b8ebf231fad876f98e7645cba8635
    
   real           beta0[N_CELLS]; // intercepts, not modeled with covariance. HMM Should they be?
   row_vector[dim] beta[N_CELLS];
   
   real<lower=0> residual_variance;  // half cauchy
 
}
    
model {    
<<<<<<< HEAD
    real logp_cov_beta[M]; // for saving log probs
    real logp_noise_neuron; // for (temporarily) saving log prob of responses given neuron is noise
    real logp_beta_neuron; // for (temporarily) saving log prob of responses given a neuron's betas
    int r; //the current response being processed, used for while loop
=======
    real logps[M]; // for saving log probs
    real logpsnoise[2]; // noise or not?
>>>>>>> c5eb4f20bc0b8ebf231fad876f98e7645cba8635
    
    // P(residual_variance)
    residual_variance ~ cauchy(0, 25);
    
<<<<<<< HEAD
    // P(mixture_weights)
    mixture_weights ~ dirichlet(alpha_cov_mix);
    
    // P(noise_weight)
    noise_weight ~ beta(alpha_noise_mix[1],alpha_noise_mix[2]);
=======
    // P(w)
    mixture_weights ~ dirichlet(alpha);
>>>>>>> c5eb4f20bc0b8ebf231fad876f98e7645cba8635
    
    // P(cov)
    for (m in 1:M) {
        covs[m] ~ inv_wishart(dim, identity);
    }
    
    // P(beta_c | covs) -- marginalize over covs for computing the betas
    for (n in 1:N_CELLS) {
        
        // put a uniform prior, and penalize below
        beta[n]  ~ uniform(-100,100);
        beta0[n] ~ normal(0,100);
        
        // but then penalize by summing over covariances
        // Directly compute the marginalized log likelihood, since stan can't handle discrete variables
        // see: http://www.michaelchughes.com/blog/2012/09/review-of-stan-off-the-shelf-hamiltonian-mcmc/
        for (m in 1:M) {
            // prob of the betas
<<<<<<< HEAD
           logp_cov_beta[m] <- log(mixture_weights[m]) + multi_normal_log(beta[n], zeros, covs[m]);
        }
        increment_log_prob(log_sum_exp(logp_cov_beta));
=======
           logps[m] <- log(mixture_weights[m]) + multi_normal_log(beta[n], zeros, covs[m]);
        }
        increment_log_prob(log_sum_exp(logps));
>>>>>>> c5eb4f20bc0b8ebf231fad876f98e7645cba8635
    }
    
    // P( y_n | beta_c) -- given those betas, how likely is the data?
    // well this comes from marginalizing over whether or not each is noise
<<<<<<< HEAD
    r <- 1;
    for(c in 1:N_CELLS) {
        logp_noise_neuron <- 0;
        logp_beta_neuron <- 0;
        while(cell[r] == c) {
            //The log prob of y, given that it is just noise. NOTE: Responses should be normalized, so using 0 should be fine
            logp_noise_neuron <- logp_noise_neuron + normal_log(y,0,residual_variance); 
        
            //The log prob of y, given betas
            logp_beta_neuron <- logp_beta_neuron + normal_log(y,beta0[cell[r]] + beta[cell[r],1]*x1[r] + beta[cell[r],2]*x2[r], residual_variance); 
        
            //y[r] ~ normal( beta0[cell[r]] + beta[cell[r],1]*x1[r] + beta[cell[r],2]*x2[r], residual_variance);
            r <- r + 1;
        }
        
        //increment log prob by the prob of this neuron's responses
        increment_log_prob(log_sum_exp(log(noise_weight) + logp_noise_neuron, log(1 - noise_weight) + logp_beta_neuron));
=======
    for(r in 1:N_RESPONSES) {
        y[r] ~ normal( beta0[cell[r]] + beta[cell[r],1]*x1[r] + beta[cell[r],2]*x2[r], residual_variance);
>>>>>>> c5eb4f20bc0b8ebf231fad876f98e7645cba8635
    }
    //print("covs=", covs); 
}

