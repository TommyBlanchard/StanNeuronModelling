
# shared R files for plotting etc. 
library(car)
library(ellipse)
library(rstan)

plot_covs <- function(myfit){
  #covs[iteration,chain,param] 
  covs = extract(myfit,pars='covs', permuted='false');
  beta = extract(myfit,pars='betas', permuted='false');
  mixture_weights = extract(myfit,pars='mixture_weights', permuted='false');
  noise_weight = extract(myfit,pars='noise_weight', permuted='false');
  scale = extract(myfit, pars='scale', permuted='false');
  
  num_covs = dim(covs)[3]/4;
  num_iter = dim(covs)[1];
  n_cells = dim(beta)[3]/2;

  whichmedian <- function(x) which.min(abs(x - median(x))) 
    plot_iter = whichmedian(covs[,1,num_covs+1]) #Find the median scale for the first cov matrix - we'll plot this
    plot(beta[plot_iter,1,1:n_cells],beta[plot_iter,1,n_cells+1:n_cells],xlab="beta for variable 1",ylab="beta for variable 2",cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
    adjusted_weight = mixture_weights[plot_iter,1,]*(1-.5*noise_weight[plot_iter]);
    for(c in 1:num_covs){
      plot_scale = scale[plot_iter,1,seq(c,length(scale[plot_iter,1,]), by=num_covs)];
      corr = matrix(covs[plot_iter,1,seq(c,length(covs[plot_iter,1,]),by=num_covs)],nrow=2);
      mycol <- col2rgb(palette()[1+c])/255      
      
      plot_bvnorm_ellipse(corr,plot_scale,adjusted_weight[c],mycol,.25)
      plot_bvnorm_ellipse(corr,plot_scale,adjusted_weight[c],mycol,.5)
      plot_bvnorm_ellipse(corr,plot_scale,adjusted_weight[c],mycol,.75)
      plot_bvnorm_ellipse(corr,plot_scale,adjusted_weight[c],mycol,.95)
      points(beta[plot_iter,1,1:n_cells],beta[plot_iter,1,n_cells+1:n_cells])
      abline(0,0,0,0)
    }
}

plot_bvnorm_ellipse <-function(corr,scale,weight,col,level){
  ell = ellipse(corr, scale,level=level)
  polygon(ell, col=rgb(col[1],col[2],col[3],weight*(1-level)))
}

plot_beta_dif <- function(myfit,given_betas){
  #Takes a stan model (regression model) and plots the difference
  #with the real betas (or supplied betas)
  #Eventually, use this to compare betas to betas obtained through
  #regular linear regression
  extracted_beta = extract(myfit,pars='betas', permuted='false');
  num_iter = dim(extracted_beta)[1];
  n_cells = dim(extracted_beta)[3]/2;
#   model_beta = matrix(extracted_beta[num_iter,1,],nrow=n_cells);
#   beta_dif = model_beta - given_betas;

  model_means = apply( extracted_beta, 3, mean) # compute the mean over samples
  plot(given_betas, model_means)
  abline(0,1)
}

hist_cov_angle <- function(myfit,ind1,ind2){
  #requires the fit, and the covariance matrix indices
  covs = extract(myfit,pars='covs', permuted='false');
  
  num_iter = dim(covs)[1];
  
  theta = numeric(length = num_iter);
  
  for(i in 1:num_iter){
    # Get the eigenvectors for the covariance matrices
      cov1 = matrix(covs[i,1,seq(ind1,length(covs[i,1,]),by=2)],nrow=2);
      eig1 = eigen(cov1)$vectors[,1];
      cov2 = matrix(covs[i,1,seq(ind2,length(covs[i,1,]),by=2)],nrow=2);
      eig2 = eigen(cov2)$vectors[,1];
      
      #Calculate angle (gives angle between -180 and +180)
      theta[i] = (atan2(cov1[2],cov1[1]) - atan2(cov2[2],cov2[1])) * (180/pi);
      
 }
  hist(theta)
  
  sorted_theta = theta[order(theta)];
  
  #print intervals and stuff
}

hist_mix_weight <- function(myfit){
  mix_weight = extract(myfit,pars='mixture_weights', permuted='false');
  #hist(mix_weight[mix_weight > .9],breaks=seq(from=.9,to=1,by=.005), col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  list_histo=hist(mix_weight,breaks=seq(from=0,to=1,by=.05), col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  x = mix_weight;
  mx = mean(x);
  n = length(x);
  n97 = round(n*.975)
  n2 = round(n*.025)
  high_int= sort(x, TRUE)[n97];
  low_int = sort(x, TRUE)[n2];
  lines( c(high_int,high_int), c(0,max(list_histo$counts)/2), col = "green", lwd = 5)
  text(high_int, max(list_histo$counts)/2 + 10 , toString(round(high_int, 2)), font=2)
  lines( c(low_int,low_int), c(0,max(list_histo$counts)/2), col = "green", lwd = 5)
  text(low_int, max(list_histo$counts/2) + 10 , toString(round(low_int, 2)), font=2)
}

hist_corr <- function(myfit){
  corr_mat = extract(myfit,pars='covs', permuted='false');
  correlation = corr_mat[,,2]
  mx = mean(correlation);
  n = dim(correlation)[1]*dim(correlation)[2];
  n97 = round(n*.975)
  n2 = round(n*.025)
  high_int= sort(correlation, TRUE)[n97];
  low_int = sort(correlation, TRUE)[n2];
  list_histo = hist(correlation,breaks=seq(from=-1,to=1,by=.05), col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  lines( c(mx,mx), c(0,max(list_histo$counts)), col = "red", lwd = 5)
  text(mx, max(list_histo$counts + 10) , paste("mean =", toString(round(mx, 2))), font=2)
  lines( c(high_int,high_int), c(0,max(list_histo$counts)/2), col = "green", lwd = 5)
  text(high_int, max(list_histo$counts)/2 + 10 , toString(round(high_int, 2)), font=2)
  lines( c(low_int,low_int), c(0,max(list_histo$counts)/2), col = "green", lwd = 5)
  text(low_int, max(list_histo$counts/2) + 10 , toString(round(low_int, 2)), font=2)
}

hist_scale <- function(myfit){
  scale = extract(myfit,pars='scale', permuted='false');
  xscale = scale[,,1]
  hist(xscale, col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  yscale = scale[,,2]
  hist(yscale, col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
}

hist_scale_ratio <- function(myfit){
  scale = extract(myfit,pars='scale', permuted='false');
  xscale = scale[,,1]
  dim(xscale) <- NULL
  yscale = scale[,,2]
  dim(yscale) <- NULL
  scaledif = xscale - yscale;
  scale_difference_ratio <- NULL;
  for(i in 1:length(xscale)){
    scale_difference_ratio[i] = scaledif[i]/(max(c(xscale[i],yscale[i])))
  }
  x = scale_difference_ratio;
  mx = mean(x);
  n = length(x);
  n97 = round(n*.975)
  n2 = round(n*.025)
  high_int= sort(x, TRUE)[n97];
  low_int = sort(x, TRUE)[n2];
  list_histo = hist(scale_difference_ratio, main = "Histogram of scale ratio", xlab="scale ratio",breaks=seq(from=-1,to=1,by=.05), col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  lines( c(mx,mx), c(0,max(list_histo$counts)), col = "red", lwd = 5)
  text(mx, max(list_histo$counts + 10) , paste("mean =", toString(round(mx, 2))), font=2)
  lines( c(high_int,high_int), c(0,max(list_histo$counts)/2), col = "green", lwd = 5)
  text(high_int, max(list_histo$counts)/2 + 10 , toString(round(high_int, 2)), font=2)
  lines( c(low_int,low_int), c(0,max(list_histo$counts)/2), col = "green", lwd = 5)
  text(low_int, max(list_histo$counts/2) + 10 , toString(round(low_int, 2)), font=2)
}

hist_noise_weight <- function(myfit){
  noise_weight = extract(myfit,pars='noise_weight', permuted='false');
  x = noise_weight;
  mx = mean(x);
  n = dim(x)[1]*dim(x)[2];
  n97 = round(n*.975)
  n2 = round(n*.025)
  high_int= sort(x, TRUE)[n97];
  low_int = sort(x, TRUE)[n2];
  list_histo = hist(noise_weight,breaks=seq(from=0,to=1,by=.05), col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  lines( c(mx,mx), c(0,max(list_histo$counts)), col = "red", lwd = 5)
  text(mx, max(list_histo$counts + 10) , paste("mean =", toString(round(mx, 2))), font=2)
  lines( c(high_int,high_int), c(0,max(list_histo$counts)/2), col = "green", lwd = 5)
  text(high_int, max(list_histo$counts)/2 + 10 , toString(round(high_int, 2)), font=2)
  lines( c(low_int,low_int), c(0,max(list_histo$counts)/2), col = "green", lwd = 5)
  text(low_int, max(list_histo$counts/2) + 10 , toString(round(low_int, 2)), font=2)
}

vis_single_cov <- function(myfit,name){
  pdf(paste(name, '_figs.pdf',sep=''))
  hist_corr(myfit)
  hist_scale(myfit)
  hist_scale_ratio(myfit)
  hist_noise_weight(myfit)
  plot_covs(myfit)
  dev.off()
}

vis_multi_cov <- function(myfit,name){
  pdf(paste(name, '_mix_hist.pdf',sep=''))
  hist_mix_weight(myfit)
  dev.off()
  pdf(paste(name, '_scatter.pdf',sep=''))
  plot_covs(myfit)
  dev.off()
}