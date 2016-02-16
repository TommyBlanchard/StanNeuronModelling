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
    matrix[2,2] identity; //Identity for scaling matrix
    
    zeros <- rep_vector(0, dim);
   
    identity <- diag_matrix(rep_vector(1.0,2));    
}

parameters {
   
   simplex[M] mixture_weights;
   
   real<lower=0> xaligned_sd; 
   real<lower=0> yaligned_sd; 
   
   vector<lower=0>[dim] free_scale; // Prior for scale (cauchy)
   corr_matrix[dim]     free_cov;
   
}
    
transformed parameters {
    vector<lower=0>[dim] xaligned_cov; 
    vector<lower=0>[dim] yaligned_cov;

    xaligned_cov[1] <- xaligned_sd;
    xaligned_cov[2] <- 0;

    yaligned_cov[1] <- yaligned_sd;
    yaligned_cov[2] <- 0;
}
    
model {    
    real logps[M]; // for saving log probs from each model
    
    mixture_weights ~ dirichlet(alpha_cov_mix);
    
    xaligned_sd ~ cauchy(0,1);
    yaligned_sd ~ cauchy(0,1);
    
    free_scale  ~ cauchy(0,1);
    free_cov    ~ lkj_corr(2);
        
    for (n in 1:N_CELLS) {
        
        // compute the probability of the betas under each component in logps
        
        // noise
        logps[1] <- multi_normal_log(X[n], zeros, sigma[n]); 
        
        // x-aligned
        logps[2] <- log(mixture_weights[2]) + multi_normal_log(X[n], zeros, diag_matrix(xaligned_cov) + sigma[n]);
        
        // y-aligned
        logps[3] <-  log(mixture_weights[2]) + multi_normal_log(X[n], zeros, diag_matrix(yaligned_cov) + sigma[n]);
        
        // free
        logps[4] <- log(mixture_weights[4]) + multi_normal_log(X[n], zeros, diag_matrix(free_scale)*free_cov*diag_matrix(free_scale) + sigma[n]);
        
        // increment by sum of probabilities, marginalizing over generating models
        increment_log_prob(log_sum_exp(logps));
    }    
}


