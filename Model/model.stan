data {
    int<lower=0> N_CELLS; // Sample size
    int<lower=1> dim; // Number of dimensions
    int<lower=1> M; // number of mixture components -- 4 here
    row_vector[dim] X[N_CELLS]; // Coefficient value for each sample on each dimension
    cov_matrix[dim] sigma[N_CELLS]; //
     
    vector[M] alpha_cov_mix; // dirichlet prior
}

transformed data {
    vector[dim] zeros; //For now, means are zero    
    zeros <- rep_vector(0, dim);
}

parameters {
   
   simplex[M] mix_weights;
   
   real<lower=0> xaligned_sd; 
   real<lower=0> yaligned_sd; 
   real<lower=0, upper=1> axes_mix_weight;
   
   real<lower=0, upper=1> noise_weight;
   
   vector<lower=0>[dim] free_scale; // Prior for scale (cauchy)
   corr_matrix[dim] free_cor;   
}
    
transformed parameters {
    vector[dim] xaligned_cov; 
    vector[dim] yaligned_cov;
    cov_matrix[dim] free_cov;

    xaligned_cov[1] <- xaligned_sd;
    xaligned_cov[2] <- 0;

    yaligned_cov[1] <- 0;
    yaligned_cov[2] <- yaligned_sd;
    
    free_cov <- diag_matrix(free_scale)*free_cor*diag_matrix(free_scale);
}
    
model {    
    real logpsignal[M]; // for saving log probs from each signal model
    real logpnoise; // for saving log prob of noise component
    real logpaxes[2]; //for saving log prob for x & y subcomponents of the axes component
    
    mix_weights ~ dirichlet(alpha_cov_mix);
    
    noise_weight ~ beta(1,1);
    axes_mix_weight ~ beta(1,1);
    
    xaligned_sd ~ cauchy(0,1);
    yaligned_sd ~ cauchy(0,1);
    
    free_scale  ~ cauchy(0,1);
    free_cor    ~ lkj_corr(1);
        
    for (n in 1:N_CELLS) {
        
        // compute the probability of the betas under each component in logps
        
        // noise
        logpnoise <- multi_normal_log(X[n], zeros, sigma[n]); 
        
        // x-aligned
        logpaxes[1] <- log(axes_mix_weight) + multi_normal_log(X[n], zeros, diag_matrix(xaligned_cov) + sigma[n]);
        
        // y-aligned
        logpaxes[2] <-  log(1 - axes_mix_weight) + multi_normal_log(X[n], zeros, diag_matrix(yaligned_cov) + sigma[n]);
        
        //xy-aligned
        logpsignal[1] <- log(mix_weights[1]) + log_sum_exp(logpaxes);
        
        // free
        logpsignal[2] <- log(mix_weights[2]) + multi_normal_log(X[n], zeros, free_cov + sigma[n]);
        
        // increment by sum of probabilities, marginalizing over generating models
        increment_log_prob(log_sum_exp(log(1 - noise_weight) + log_sum_exp(logpsignal), log(noise_weight) + logpnoise));
    }    
}


