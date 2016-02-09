data {
    int<lower=0> N_CELLS; // Sample size
    int<lower=1> dim; // Number of dimensions
    row_vector[dim] x[N_CELLS]; // Coefficient value for each sample on each dimension
    cov_matrix[dim] sigma[N_CELLS]; //covariance matrix for coefficients
    
    //submodels: 1 free, 2 free, 1 circle (no correlation), 'plus' (along axes)
    vector[2] alpha_model_mix; // dirichlet prior//vector[4] alpha_model_mix; // dirichlet prior
    vector[2] alpha_noise_mix; //beta prior
}

transformed data {
    vector[dim] zeros; //For now, means are zero
    corr_matrix[dim] no_corr;
    
    zeros <- rep_vector(0, dim);
    no_corr <- diag_matrix(rep_vector(1,dim));
}

parameters {
    //1 free submodel parameters
    corr_matrix[dim] cor1free;
    vector<lower=0>[dim] scale1free;
    row_vector[dim] betas1free[N_CELLS];
    
    //2 covs submodel parameters
    corr_matrix[dim] cor2free[2];
    vector<lower=0>[dim] scale2free[2];
    simplex[2] mix2free;
    row_vector[dim] betas2free[N_CELLS];
    
    // no corr submodel parameters
    //vector<lower=0>[dim] scale1circ;
    //row_vector[dim] betas1circ[N_CELLS];
    
    // axes submodel parameters
    //vector<lower=0>[dim] scale2axes;
    //simplex[2] mix2axes;
    //row_vector[dim] betas2axes[N_CELLS];
       
    //Shared noise weight component and beta estimates
    real<lower=0, upper=1> noise_weight;

   
   //model mixture weights
   simplex[2] model_mixture_weights;//simplex[4] model_mixture_weights;
}

//transformed  parameters{
//    row_vector[dim] betas2axes_dist1[N_CELLS];
//    row_vector[dim] betas2axes_dist2[N_CELLS];
//    
//    for (n in 1:N_CELLS) {
//        betas2axes_dist1[n,1] <- betas2axes[n,1];
//        betas2axes_dist1[n,2] <- 0;
//        betas2axes_dist2[n,2] <- betas2axes[n,2];
//        betas2axes_dist2[n,1] <- 0;
//    }
//}
    
model {    
    real logp2free_parts[2];
    //real logp2axes_parts[2];
    
    real logp_submodels[2];//[4];
    
    real logp_betas[2];//[4];
    real logp_data[2];//[4];
    
    real logp_noise; //for saving log prob of obs coming from noise
    real logp_signal; //for saving log prob of obs coming from signal
        
    // P(model_mixture_weights)
    model_mixture_weights ~ dirichlet(alpha_model_mix);
    
    // P(noise_weight)
    noise_weight ~ beta(alpha_noise_mix[1],alpha_noise_mix[2]);
    
    // PRIORS FOR ALL OF THE SUBMODELS 
    
    // P(corr, scaling free) 1 free model
    scale1free ~ cauchy(0, 2.5);
    cor1free ~ lkj_corr(1); //lkj(1) gives equal prob to positive, negative, and 0 correlation

    // P(corr, scaling free, mixture) 2 free model
    for (m in 1:2){
        scale2free[m] ~ cauchy(0, 2.5);
        cor2free[m] ~ lkj_corr(1); //lkj(1) gives equal prob to positive, negative, and 0 correlation
    }
    mix2free ~ beta(1,1); //flat prior for mixture weights
    
    // P(scale) no correlation model
    //scale1circ ~ cauchy(0, 2.5);
    
    // P(scale, mixture) axes model
    //scale2axes ~ cauchy(0, 2.5);    
    //mix2axes ~ beta(1,1); //flat prior for mixture weights
    
    //P(DATA | FULL MODEL)
    for (n in 1:N_CELLS) {
        //P(DATA | EACH SUBMODEL)
        
        //p(betas | 1free params)
        logp_betas[1] <- multi_normal_log(betas1free[n], zeros, diag_matrix(scale1free)*cor1free*diag_matrix(scale1free));
        //p(data | betas)
        logp_data[1] <- multi_normal_log(x[n], betas1free[n], sigma[n]);
        logp_submodels[1] <- log(model_mixture_weights[1]) + log_sum_exp(logp_betas[1], logp_data[1]);
        
        //p(betas | 2 free dists)
        for (m in 1:2) {
            logp2free_parts[m] <- mix2free[m] + multi_normal_log(betas2free[n], zeros, diag_matrix(scale2free[m])*cor2free[m]*diag_matrix(scale2free[m]));
        }
        logp_betas[2] <- log_sum_exp(logp2free_parts);
        logp_data[2] <- multi_normal_log(x[n], betas2free[n], sigma[n]);
        logp_submodels[2] <- log(model_mixture_weights[2]) + log_sum_exp(logp_betas[2], logp_data[2]);
        
        //p(betas | big circle)
        //logp_betas[3] <- multi_normal_log(betas1circ[n], zeros, diag_matrix(scale1circ)*no_corr*diag_matrix(scale1circ));        
        //logp_data[3] <-  multi_normal_log(x[n], betas1circ[n], sigma[n]);
        //logp_submodels[3] <- log(model_mixture_weights[3]) + log_sum_exp(logp_betas[3], logp_data[3]);
        
        //p(betas | axes)
        //logp_betas[4] <- multi_normal_log(betas2axes[n], zeros, diag_matrix(scale2axes)*no_corr*diag_matrix(scale2axes));
        //logp2axes_parts[1] <- mix2axes[1] + multi_normal_log(x[n], betas2axes_dist1[n], sigma[n]);     
        //logp2axes_parts[2] <- mix2axes[2] + multi_normal_log(x[n], betas2axes_dist2[n], sigma[n]);
        //logp_data[4] <- log_sum_exp(logp2axes_parts);
        //logp_submodels[4] <- log(model_mixture_weights[4]) + log_sum_exp(logp_betas[4], logp_data[4]);
        
        //p(observation | this cell is signal)
        logp_signal <- log_sum_exp(logp_submodels);
        
        //p(observation | this cell is noise)
        logp_noise <- multi_normal_log(x[n], zeros, sigma[n]);        
        
        //increment prob by p(obs | noise, signal)
        increment_log_prob(log_sum_exp(log(1.0 - noise_weight) + logp_signal, log(noise_weight) + logp_noise));
    }    
}