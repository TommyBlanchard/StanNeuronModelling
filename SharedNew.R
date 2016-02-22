
# shared R files for plotting etc. 
library(car)
library(ellipse)
library(rstan)

plot_rand_iter <- function(myfit, data){
  pars <- extract(myfit, permuted='true');
  num_iter <- dim(pars$free_cov)[1];
  plot_iter <- sample(1:num_iter, 1)
  plot_covs(myfit,data,plot_iter)
}

plot_covs <- function(myfit, data, plot_iter){
  #covs[iteration,chain,param] 
  pars = extract(myfit, permuted='true')
  mean_sigma = colMeans(data$sigma);
  
  plot(data$X, pch=20,ylim=c(min(data$X),max(data$X)),xlim=c(min(data$X),max(data$X)),xlab="beta for variable 1",ylab="beta for variable 2",cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  
  cov_to_plot = diag(pars$xaligned_cov[plot_iter,]) + mean_sigma;
  mycol <- col2rgb(palette()[3])/255;
  plot_ellipse_dens(cov_to_plot,pars$mix_weights[plot_iter,1],mycol)
  
  cov_to_plot = diag(pars$yaligned_cov[plot_iter,]) + mean_sigma;
  mycol <- col2rgb(palette()[3])/255;
  plot_ellipse_dens(cov_to_plot,pars$mix_weights[plot_iter,1],mycol)
  
  cov_to_plot = pars$free_cov[plot_iter,,] + mean_sigma;
  mycol <- col2rgb(palette()[4])/255;
  plot_ellipse_dens(cov_to_plot,pars$mix_weights[plot_iter,2],mycol)
  
  cov_to_plot = mean_sigma;
  mycol <- col2rgb(palette()[2])/255;
  plot_ellipse_dens(cov_to_plot,pars$noise_weight[plot_iter],mycol)
  
  points(data$X, pch=20)
  
  abline(0,0,0,0)
}

plot_ellipse_dens <-function(cov,weight,mycol){
  plot_bvnorm_ellipse(cov,weight,mycol,.5)
  plot_bvnorm_ellipse(cov,weight,mycol,.95)
}

plot_bvnorm_ellipse <-function(cov,weight,col,level){
  ell = ellipse(cov, level=level)
  polygon(ell, border=rgb(col[1],col[2],col[3]),col=rgb(col[1],col[2],col[3],weight*(1-level/2))) #level/2 here just to make sure the higher levels still get some color!
}

hist_mix_weight <- function(myfit){
  #plots the historgram of the mix weight between the free and axes submodels
  pars = extract(myfit, permuted='true')
  mix_weight = extract(myfit,pars='mix_weights', permuted='false');
  #hist(mix_weight[mix_weight > .9],breaks=seq(from=.9,to=1,by=.005), col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  list_histo=hist(mix_weight[mix_weight < .5],breaks=seq(from=0,to=0.5,by=.01), col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  x = mix_weight[mix_weight < .5];
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
  pars = extract(myfit, permuted='true')
  
  correlation = pars$free_cor[,1,2]
  mx = mean(correlation);
  n = length(correlation)[1];
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

hist_free_scale <- function(myfit){
  pars = extract(myfit, permuted='true')
  
  xscale = pars$free_scale[,1]
  hist(xscale, col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  yscale = pars$free_scale[,2]
  hist(yscale, col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
}

hist_axes_scale <- function(myfit){
  pars = extract(myfit, permuted='true')
  
  xscale = pars$xaligned_sd
  hist(xscale, col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  yscale = pars$xaligned_sd
  hist(yscale, col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
}

free_scale_ratio <- function(myfit){
  pars = extract(myfit, permuted='true')
  xscale = pars$free_scale[,1]
  dim(xscale) <- NULL
  yscale = pars$free_scale[,2]
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
  list_histo = hist(scale_difference_ratio, main = "Histogram of free scale ratio", xlab="scale ratio",breaks=seq(from=-1,to=1,by=.05), col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  lines( c(mx,mx), c(0,max(list_histo$counts) - 10), col = "red", lwd = 5)
  text(mx, max(list_histo$counts) , paste("mean =", toString(round(mx, 2))), font=2)
  lines( c(high_int,high_int), c(0,max(list_histo$counts)/2), col = "green", lwd = 5)
  text(high_int, max(list_histo$counts)/2 + 10 , toString(round(high_int, 2)), font=2)
  lines( c(low_int,low_int), c(0,max(list_histo$counts)/2), col = "green", lwd = 5)
  text(low_int, max(list_histo$counts/2) + 10 , toString(round(low_int, 2)), font=2)
}

hist_noise_weight <- function(myfit){
  pars = extract(myfit, permuted='true')
  noise_weight = pars$noise_weight;
  mx = mean(noise_weight);
  n = length(noise_weight);
  n97 = round(n*.975)
  n2 = round(n*.025)
  high_int= sort(noise_weight, TRUE)[n97];
  low_int = sort(noise_weight, TRUE)[n2];
  list_histo = hist(noise_weight,breaks=seq(from=0,to=1,by=.05), col='blue',cex.lab=2, cex.main=2, cex.axis=2, font.lab=2, font.main=2, font.axis=2)
  lines( c(mx,mx), c(0,max(list_histo$counts) - 10), col = "red", lwd = 5)
  text(mx, max(list_histo$counts) , paste("mean =", toString(round(mx, 2))), font=2)
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