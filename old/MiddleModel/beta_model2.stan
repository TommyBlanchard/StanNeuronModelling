data {
    int<lower=0> N_CELLS; // Sample size
    int<lower=1> dim; // Number of dimensions
    int<lower=1> M; // number of mixture components -- 4 here
    row_vector[dim] X[N_CELLS]; // Coefficient value for each sample on each dimension
    
    vector[2] xWishartPrior; // the priors for x and y wisharts
    vector[2] yWishartPrior;
        
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
   
   vector<lower=0>[dim] noise_scale; // Prior for scale (cauchy)
   vector<lower=0>[dim] noise_sd; // Prior for scale (cauchy)
   
   vector<lower=0>[dim] xaligned_scale; // Prior for scale (cauchy)
   cov_matrix[dim]      xaligned_cov;
   
   vector<lower=0>[dim] yaligned_scale; // Prior for scale (cauchy)
   cov_matrix[dim]      yaligned_cov;
   
   vector<lower=0>[dim] free_scale; // Prior for scale (cauchy)
   corr_matrix[dim]     free_cov;
   
}
    
model {    
    real logps[M]; // for saving log probs from each model
    
    mixture_weights ~ dirichlet(alpha_cov_mix);
    
    noise_scale    ~ cauchy(0,1); // TODO: MAYBE THE SCALE SHOULD BE DRAWN FROM A WISHART?
	noise_sd       ~ cauchy(0,1);
    
    xaligned_scale ~ cauchy(0,1);
	xaligned_cov   ~ inv_wishart(dim, diag_matrix(xWishartPrior));
	
    yaligned_scale ~ cauchy(0,1);
	yaligned_cov   ~ inv_wishart(dim, diag_matrix(yWishartPrior));
	
	free_scale     ~ cauchy(0,1);
	free_cov       ~ inv_wishart(dim, identity);
        
    for (n in 1:N_CELLS) {
        
		// probability of the betas under each component
		logps[1] <- log(mixture_weights[1]) + multi_normal_log(X[n], zeros, diag_matrix(noise_scale)*diag_matrix(noise_sd)*diag_matrix(noise_scale));     
        
		logps[2] <- log(mixture_weights[2]) + multi_normal_log(X[n], zeros, diag_matrix(xaligned_scale)*xaligned_cov*diag_matrix(xaligned_scale));
        
		logps[3] <- log(mixture_weights[3]) + multi_normal_log(X[n], zeros, diag_matrix(yaligned_scale)*yaligned_cov*diag_matrix(yaligned_scale));
        
		logps[4] <- log(mixture_weights[4]) + multi_normal_log(X[n], zeros, diag_matrix(free_scale)*free_cov*diag_matrix(free_scale));
		
		increment_log_prob(log_sum_exp(logps));
    }    
}
