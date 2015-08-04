
# shared R files for plotting etc. 
library(car)

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
  scale = extract(myfit, pars='scale', permuted='false');
  
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
    for(c in 1:num_covs){
      unscaled_cov = matrix(covs[plot_iter,1,seq(c,length(covs[plot_iter,1,]),by=2)],nrow=2);
      s = diag(c(scale[plot_iter,1,1], scale[plot_iter,1,2]))
      cur_cov = s %*% unscaled_cov %*% s
      
      e <- eigen(cur_cov)
      
    # From http://stats.stackexchange.com/questions/9898/how-to-plot-an-ellipse-from-eigenvalues-and-eigenvectors-in-r
    ctr <- c(0,0)
    angles <- seq(0, 2*pi, length.out=200) 
      eigVal  <- eigen(cur_cov)$values
    eigVec  <- eigen(cur_cov)$vectors
    eigScl  <- eigVec %*% diag(sqrt(eigVal))  # scale eigenvectors to length = square-root
    xMat    <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
    yMat    <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])
    ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles)) # normal ellipse
    ellRot  <- eigVec %*% t(ellBase)  # rotated ellipse
    mycol <- palette()[1+c]
    lines((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, type="l", lwd=2, col=mycol)
    matlines(xMat, yMat, lty=1, lwd=2, col=mycol)
    points(ctr[1], ctr[2], pch=4, col=mycol, lwd=3)
      
#       ellipse(c(0, 0), shape=cur_cov, radius=1, col=palette()[1+c], fill=TRUE, fill.alpha = adjusted_weight[c])
    }
  }
}

plot_beta_dif <- function(myfit,given_betas){
  #Takes a stan model (regression model) and plots the difference
  #with the real betas (or supplied betas)
  #Eventually, use this to compare betas to betas obtained through
  #regular linear regression
  extracted_beta = extract(myfit,pars='beta', permuted='false');
  num_iter = dim(extracted_beta)[1];
  n_cells = dim(extracted_beta)[3]/2;
#   model_beta = matrix(extracted_beta[num_iter,1,],nrow=n_cells);
#   beta_dif = model_beta - given_betas;

  model_means = apply( extracted_beta, 3, mean) # compute the mean over samples
  plot(given_betas, model_means)
  abline(0,1)
}