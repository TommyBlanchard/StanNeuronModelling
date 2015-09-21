// List of matrices in R: https://groups.google.com/forum/#!topic/stan-users/-6s9PvVsyug

data {
    int<lower=0> N_CELLS; // Sample size
    int<lower=1> dim; // Number of dimensions
    int<lower=1> M; // number of mixture components 
    row_vector[dim] x[N_CELLS]; // Coefficient value for each sample on each dimension
    vector[dim] sigma[N_CELLS]; //standard error for coefficients
    
    vector[M] alpha_cov_mix; // dirichlet prior
    vector[2] alpha_noise_mix; //beta prior
}

transformed data {
    vector[dim] zeros; //For now, means are zero
    matrix[dim,dim] identity; //Identity for scaling matrix
    
    zeros <- rep_vector(0, dim);
    identity <- diag_matrix(rep_vector(1.0,dim)); 
}

parameters {
   corr_matrix[dim] covs[M]; //If we keep the lkj version, rename this variable - no longer a covariance matrix!
   
   simplex[M] mixture_weights;
   real<lower=0, upper=1> noise_weight;
   
   row_vector[dim] beta[N_CELLS];
      
   vector<lower=0>[dim] scale[M]; // Prior for scale (cauchy)
}
    
model {    
    real logps[M]; // for saving log probs
    real logpnoise; // for saving log prob of beta coming from noise mixture
    
    real logp_noise; //for saving log prob of obs coming from noise
    real logp_signal; //for saving log prob of obs coming from signal
    
    // P(mixture_weights)
    mixture_weights ~ dirichlet(alpha_cov_mix);
    
    // P(noise_weight)
    noise_weight ~ beta(alpha_noise_mix[1],alpha_noise_mix[2]);
    
    // P(cov), P(scaling factors)
    for (m in 1:M) {
        scale[m] ~ cauchy(0, 1.5);
        covs[m] ~ lkj_corr(1); //lkj(1) gives equal prob to positive, negative, and 0 correlation
    }
        
    for (n in 1:N_CELLS) {
        //p(betas | cov)
        for (m in 1:M) {
            logps[m] <- log(mixture_weights[m]) + multi_normal_log(beta[n], zeros, diag_matrix(scale[m])*covs[m]*diag_matrix(scale[m]));
        }
        //p(observation | this cell is signal)
        logp_signal <- multi_normal_log(x[n], beta[n], diag_matrix(sigma[n])) + log_sum_exp(logps);
        
        //p(observation | this cell is noise)
        logp_noise <- multi_normal_log(x[n], zeros, diag_matrix(sigma[n]));
        
        //increment prob by p(obs | noise or signal)
        increment_log_prob(log_sum_exp(log(1.0 - noise_weight) + logp_signal, log(noise_weight) + logp_noise));
    }    
}