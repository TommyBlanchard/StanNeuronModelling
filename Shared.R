
# shared R files for plotting etc. 


plot_covs <- function(myfit,num_plots){
  #Plots the betas and draws ellipses to represent the covariance matrices
  #myfit is the fit given by stan
  #num_plots is the number of plots you want to see: a value of  2 will
  #show betas & covariances matrices for halfway through the iterations 
  #and the last iteration.
  
  #TODO: currently will only work with the regression model.
  # currently plots for just 1 chain, ignores others
  
  #covs[iteration,chain,param] 
  covs = extract(myfit,pars='covs', permuted='false');
  beta = extract(myfit,pars='beta', permuted='false');
  mixture_weights = extract(myfit,pars='mixture_weights', permuted='false');
  noise_weight = extract(myfit,pars='noise_weight', permuted='false');
  
  num_covs = dim(covs)[3]/4;
  num_iter = dim(covs)[1];
  n_cells = dim(beta)[3]/2;
  
  iter = floor(num_iter/num_plots);
  
  plot_rows = sqrt(num_plots);
  plot_cols = ceiling(num_plots/plot_rows);
  par(mfrow=c(plot_rows,plot_cols))
  for(i in 1:num_plots){
    print(i)
    plot_iter = iter*i;
    plot(beta[plot_iter,1,1:n_cells],beta[plot_iter,1,n_cells+1:n_cells])
    abline(h=0,v=0)
    adjusted_weight = (1-noise_weight[plot_iter])*mixture_weights[plot_iter,1,];
    for(c in 1:numcovs){
      cur_cov = matrix(covs[plot_iter,1,seq(c,length(covs[plot_iter,1,]),by=2)],nrow=2);
      ellipse(c(0, 0), shape=cur_cov, radius=1, col=palette()[1+c], fill=TRUE, fill.alpha = adjusted_weight[c])
    }
  }
}

plot_beta_dif <- function(myfit,given_betas){
  #Takes a stan model (regression model) and plots the difference
  #with the real betas (or supplied betas)
  #Eventually, use this to compare betas to betas obtained through
  #regular linear regression
  extracted_beta = extract(myfit,pars='beta', permuted='false');
  num_iter = dim(beta)[1];
  n_cells = dim(beta)[3]/2;
  model_beta = matrix(model_beta[num_iter,1,],nrow=n_cells);
  beta_dif = model_beta - given_betas;
  plot(beta_dif)
}