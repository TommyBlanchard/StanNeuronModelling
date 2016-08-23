
# shared R files for plotting etc. 
library(car)
library(ellipse)
library(rstan)
library(coda)
library(ggplot2)
library(gridExtra)


plot_rand_iter <- function(myfit, data){
  pars <- extract(myfit, permuted='true');
  num_iter <- dim(pars$free_cov)[1];
  plot_iter <- sample(1:num_iter, 1)
  plot_covs(myfit,data,plot_iter)
}

plot_data_covs <- function(data){
  plot(data$X, main="Data covariances",  pch=20,ylim=c(-max(abs(c(data$X))),max(abs(c(data$X)))),xlim=c(-max(abs(c(data$X))),max(abs(c(data$X)))),cex.main=1,cex.axis=1.1, font.lab=2, font.main=2, font.axis=2)
  abline(0,0,0,0)
  for(i in 1:dim(data$X)[1]){
    plot_bvnorm_ellipse(data$sigma[i,,],1,col2rgb(palette()[8])/255,.95, data$X[i,])
  }
  points(data$X, pch=20)
}

plot_5free_iters<- function(myfit,data){
  #covs[iteration,chain,param] 
  pars = extract(myfit, permuted='false')
  num_iter <- dim(pars)[1]; #Number of iters per chain
  pars = extract(myfit, permuted='true')
  
  plot(data$X, main="Mixed tuning distributions",  pch=20,ylim=c(-max(abs(c(data$X))),max(abs(c(data$X)))),xlim=c(-max(abs(c(data$X))),max(abs(c(data$X)))),cex.main=1,cex.axis=1.1, font.lab=2, font.main=2, font.axis=2)
  abline(0,0,0,0)
  for( i in 1:5){
    plot_iter = round(num_iter*(i/5))
    cov_to_plot = pars$free_cov[plot_iter,,];
    mycol <- col2rgb(palette()[4])/255;
    plot_ellipse_dens(cov_to_plot,pars$mix_weights[plot_iter,2],mycol)
  }
  
  logpnoise = colMeans(pars$logpnoise)
  logpfree = colMeans(pars$logpfree)
  logpaxes = colMeans(log(exp(pars$logpyaxis) + exp(pars$logpxaxis)))
  
  plot_points(data$X,logpnoise,logpfree,logpaxes, col2rgb(palette()[2])/255) #plot noise points
  plot_points(data$X,logpfree,logpaxes,logpnoise, col2rgb(palette()[4])/255) #plot axes points
  plot_points(data$X,logpaxes,logpnoise,logpfree, col2rgb(palette()[3])/255) #plot free points
  

}
  
plot_5axes_iters<- function(myfit,data){
  #covs[iteration,chain,param] 
  pars = extract(myfit, permuted='false')
  num_iter <- dim(pars)[1]; #Number of iters per chain
  pars = extract(myfit, permuted='true')
  
  plot(data$X, main="Pure tuning distributions",  pch=20,ylim=c(-max(abs(c(data$X))),max(abs(c(data$X)))),xlim=c(-max(abs(c(data$X))),max(abs(c(data$X)))),cex.main=1,cex.axis=1.1, font.lab=2, font.main=2, font.axis=2)
  abline(0,0,0,0)
  for( i in 1:5){
    plot_iter = round(num_iter*(i/5))
    
    #Find 10 points most likely to be from x_aligned cov
    mean_sigma = get_mean_sigma(pars$logpxaxis[plot_iter,],data)
    cov_to_plot = diag(pars$xaligned_cov[plot_iter,]) + mean_sigma;
    mycol <- col2rgb(palette()[3])/255;
    plot_ellipse_dens(cov_to_plot,pars$mix_weights[plot_iter,1],mycol)
    
    mean_sigma = get_mean_sigma(pars$logpyaxis[plot_iter,],data)
    cov_to_plot = diag(pars$yaligned_cov[plot_iter,]) + mean_sigma;
    mycol <- col2rgb(palette()[3])/255;
    plot_ellipse_dens(cov_to_plot,pars$mix_weights[plot_iter,1],mycol)
  }    
  
  logpnoise = colMeans(pars$logpnoise)
  logpfree = colMeans(pars$logpfree)
  logpaxes = colMeans(log(exp(pars$logpyaxis) + exp(pars$logpxaxis)))
  
  plot_points(data$X,logpnoise,logpfree,logpaxes, col2rgb(palette()[2])/255) #plot noise points
  plot_points(data$X,logpfree,logpaxes,logpnoise, col2rgb(palette()[4])/255) #plot axes points
  plot_points(data$X,logpaxes,logpnoise,logpfree, col2rgb(palette()[3])/255) #plot free points
  
  
}
 
plot_5noise_iters<- function(myfit,data){
  #covs[iteration,chain,param] 
  pars = extract(myfit, permuted='false')
  num_iter <- dim(pars)[1]; #Number of iters per chain
  pars = extract(myfit, permuted='true')
  
  plot(data$X, main="Noise distributions",  pch=20,ylim=c(-max(abs(c(data$X))),max(abs(c(data$X)))),xlim=c(-max(abs(c(data$X))),max(abs(c(data$X)))), cex.main=1,cex.axis=1.1, font.lab=2, font.main=2, font.axis=2)
  abline(0,0,0,0)
  for( i in 1:5){
    plot_iter = round(num_iter*(i/5))
    mean_sigma = get_mean_sigma(pars$logpnoise[plot_iter,],data)
    mm_sigma = mean(mean_sigma[1,1],mean_sigma[2,2])
    cov_to_plot = diag(c(mm_sigma,mm_sigma));
    mycol <- col2rgb(palette()[2])/255;
    plot_ellipse_dens(cov_to_plot,pars$noise_weight[plot_iter],mycol)
  }
  
  logpnoise = colMeans(pars$logpnoise)
  logpfree = colMeans(pars$logpfree)
  logpaxes = colMeans(log(exp(pars$logpyaxis) + exp(pars$logpxaxis)))
  
  plot_points(data$X,logpnoise,logpfree,logpaxes, col2rgb(palette()[2])/255) #plot noise points
  plot_points(data$X,logpfree,logpaxes,logpnoise, col2rgb(palette()[4])/255) #plot axes points
  plot_points(data$X,logpaxes,logpnoise,logpfree, col2rgb(palette()[3])/255) #plot free points
}
 
plot_5cov_iters <- function(myfit,data){
  #covs[iteration,chain,param] 
  pars = extract(myfit, permuted='false')
  num_iter <- dim(pars)[1]; #Number of iters per chain
  pars = extract(myfit, permuted='true')
  
  plot(data$X, main="Full model distributions", pch=20,ylim=c(-max(abs(c(data$X))),max(abs(c(data$X)))),xlim=c(-max(abs(c(data$X))),max(abs(c(data$X)))),cex.main=1,cex.axis=1.1, font.lab=2, font.main=2, font.axis=2)
  abline(0,0,0,0)
  for( i in 1:5){
    plot_iter = round(num_iter*(i/5))
    
    #Find 10 points most likely to be from x_aligned cov
    mean_sigma = get_mean_sigma(pars$logpxaxis[plot_iter,],data)
    cov_to_plot = diag(pars$xaligned_cov[plot_iter,]) + mean_sigma;
    mycol <- col2rgb(palette()[3])/255;
    plot_ellipse_dens(cov_to_plot,pars$mix_weights[plot_iter,1],mycol)
    
    mean_sigma = get_mean_sigma(pars$logpyaxis[plot_iter,],data)
    cov_to_plot = diag(pars$yaligned_cov[plot_iter,]) + mean_sigma;
    mycol <- col2rgb(palette()[3])/255;
    plot_ellipse_dens(cov_to_plot,pars$mix_weights[plot_iter,1],mycol)
  }    
  for( i in 1:5){
    plot_iter = round(num_iter*(i/5))
    cov_to_plot = pars$free_cov[plot_iter,,];
    mycol <- col2rgb(palette()[4])/255;
    plot_ellipse_dens(cov_to_plot,pars$mix_weights[plot_iter,2],mycol)
  }
  for( i in 1:5){
    plot_iter = round(num_iter*(i/5))
    mean_sigma = get_mean_sigma(pars$logpnoise[plot_iter,],data)
    mm_sigma = mean(mean_sigma[1,1],mean_sigma[2,2])
    cov_to_plot = diag(c(mm_sigma,mm_sigma));
    mycol <- col2rgb(palette()[2])/255;
    plot_ellipse_dens(cov_to_plot,pars$noise_weight[plot_iter],mycol)
  }

  logpnoise = colMeans(pars$logpnoise)
  logpfree = colMeans(pars$logpfree)
  logpaxes = colMeans(log(exp(pars$logpyaxis) + exp(pars$logpxaxis)))
  
  plot_points(data$X,logpnoise,logpfree,logpaxes, col2rgb(palette()[2])/255) #plot noise points
  plot_points(data$X,logpfree,logpaxes,logpnoise, col2rgb(palette()[4])/255) #plot axes points
  plot_points(data$X,logpaxes,logpnoise,logpfree, col2rgb(palette()[3])/255) #plot free points
}

get_mean_sigma <- function(logp,data){
  neurons = unlist(as.matrix(sort.int(logp, index.return=TRUE,decreasing=TRUE))[2])[1:10]
  mean_sigma = colMeans(data$sigma[neurons,,])
  return(mean_sigma)
}

plot_ellipse_dens <-function(cov,weight,mycol){
  plot_bvnorm_ellipse(cov,weight,mycol,.95, c(0,0))
}

plot_bvnorm_ellipse <-function(cov,weight,mycol,level, centre, bordercol=NULL){
  if (is.null(bordercol)){
    bordercol = mycol;
  }
  ell = ellipse(cov, centre = centre, level=level)
  polygon(ell,col=rgb(mycol[1],mycol[2],mycol[3],weight/5),border=rgb(bordercol[1],bordercol[2],bordercol[3])) 
}

plot_points <-function(X,logp1,logp2,logp3,mycol){ #plot points X in col=col if logp1 > logp2, logp3
  ind = (logp1 > logp2) & (logp1 > logp3);
  points(X[ind,],pch=20, col=rgb(mycol[1],mycol[2],mycol[3]))
}

hist_mix_weight <- function(myfit){
  #plots the historgram of the mix weight between the free and axes submodels
  pars = extract(myfit, permuted='true')
  x = pars$mix_weights[,2];
  list_histo = hist(x,breaks=seq(from=0,to=1,by=.05), plot=FALSE)
  plot_hist(x,list_histo,"Posterior mixed-tuning signal weight")
}

hist_corr <- function(myfit){
  pars = extract(myfit, permuted='true')
  x = pars$free_cor[,1,2]
  list_histo = hist(x,plot=FALSE,breaks=seq(from=-1,to=1,by=.1))
  plot_hist(x,list_histo,"Posterior correlation")
}

free_scale_ratio <- function(myfit){
  pars = extract(myfit, permuted='true')
  xscale = pars$free_scale[,1]
  yscale = pars$free_scale[,2]
  scaledif = xscale - yscale;
  scale_difference_ratio <- NULL;
  for(i in 1:length(xscale)){
    scale_difference_ratio[i] = scaledif[i]/(max(c(xscale[i],yscale[i])))
  }
  x = array(scale_difference_ratio);
  list_histo = hist(scale_difference_ratio, breaks=seq(from=-1,to=1,by=.1), plot=FALSE)
  plot_hist(x,list_histo,"Posterior mixed-tuning scale difference ratio")
}

axes_scale_ratio <- function(myfit){
  pars = extract(myfit, permuted='true')
  xscale = pars$xaligned_sd
  yscale = pars$yaligned_sd
  scaledif = xscale - yscale;
  scale_difference_ratio <- NULL;
  for(i in 1:length(xscale)){
    scale_difference_ratio[i] = scaledif[i]/(max(c(xscale[i],yscale[i])))
  }
  x = array(scale_difference_ratio);
  list_histo = hist(scale_difference_ratio, breaks=seq(from=-1,to=1,by=.1), plot=FALSE)
  plot_hist(x,list_histo,"Posterior pure-tuning scale difference ratio")
}

hist_axes_weight <- function(myfit){
  pars = extract(myfit, permuted='true')
  x = array(pars$axes_mix_weight);
  list_histo = hist(x,plot=FALSE,breaks=seq(from=0,to=1,by=.05))
  plot_hist(x,list_histo,"Posterior pure-tuning variable 1 weight")
}

hist_noise_weight <- function(myfit){
  pars = extract(myfit, permuted='true')
  x = array(pars$noise_weight);
  list_histo = hist(x,plot=FALSE,breaks=seq(from=0,to=1,by=.05))
  plot_hist(x,list_histo,"Posterior noise weight")
}

hist_free_scalex <- function(myfit){
  pars = extract(myfit, permuted='true')
  x = array(pars$free_scale[,1]);
  list_histo = hist(remove_outliers(x),plot=FALSE)
  plot_hist(x,list_histo,"Posterior mixed-tuning x scale")
}

hist_free_scaley <- function(myfit){
  pars = extract(myfit, permuted='true')
  x = array(pars$free_scale[,2]);
  list_histo = hist(remove_outliers(x),plot=FALSE)
  plot_hist(x,list_histo,"Posterior mixed-tuning y scale")
}

hist_axes_scalex <- function(myfit){
  pars = extract(myfit, permuted='true')
  x = array(pars$xaligned_sd);
  list_histo = hist(remove_outliers(x),plot=FALSE)
  plot_hist(x,list_histo,"Posterior pure-tuning x scale")
}

hist_axes_scaley <- function(myfit){
  pars = extract(myfit, permuted='true')
  x = array(pars$yaligned_sd);
  list_histo = hist(remove_outliers(x),plot=FALSE)
  plot_hist(x,list_histo,"Posterior pure-tuning y scale")
}

plot_hist <- function(x,histo,main){
  mx = median(x);
  high= HPDinterval(as.mcmc(x), prob=0.95)[1];
  low = HPDinterval(as.mcmc(x), prob=0.95)[2];
  histo$counts=histo$counts/sum(histo$counts)
  plot(histo, col='grey',cex.lab=1.5, cex.main=1.5, cex.axis=1.5, font.lab=2, font.main=2, font.axis=1.5, border="#777777",main=NULL,xlab=NULL,ylab=NULL)
  title(main=main, ylab="Probability density")
  lines( c(mx,mx), c(0,par('usr')[4]*.65), col = "red", lwd = 2)
  text(mx, par('usr')[4]*.75, paste(toString(round(mx, 2))), font=2)
  
  lines( c(low,high), c(0,0), col = "red", lwd = 2)
  text(low, par('usr')[4]*.1, paste(toString(round(low, 2))), font=2)
  text(high, par('usr')[4]*.1, paste(toString(round(high, 2))), font=2)
}

#VERY lenient outlier remover only for scale plots, not calculations of HPD intervals or medians
#Without this, the scale can sometimes go off into crazy land for a component with very little weight, ruining the hist
remove_outliers <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75))
  H <- 10 * IQR(x)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

vis_dists <- function(myfit,data,X1name, X2name, filename){
  pdf(paste(filename, '_dists.pdf',sep=''), width=4,height=7.5)
  layout(matrix(c(1,2,3,4,5,5,5,5), 4, byrow='true'), widths = rep.int(1, 2),
         heights = rep.int(1, 4), respect = TRUE)
  par(pty="s",
      mar = c(2, 1, 1, 0) + 0.1,
      oma = c(2,2,0,0) + 0.1)
  plot_data_covs(data)
  plot_5free_iters(myfit,data)
  plot_5axes_iters(myfit,data)
  plot_5noise_iters(myfit,data)
  plot_5cov_iters(myfit,data)
  title(xlab = paste(X1name, 'tuning'),
        ylab = paste(X2name, 'tuning'),
        outer = TRUE,
        line=0.1,
        font.lab=2,
        cex.lab = 1.5)
  
  dev.off()
}

vis_posteriors <- function(myfit,name){
  pdf(paste(name, '_posteriors.pdf',sep=''))
  par(mfrow=c(5,2), mar=c(2,2.5,2,1))
  hist_mix_weight(myfit)
  hist_corr(myfit)
  free_scale_ratio(myfit)
  axes_scale_ratio(myfit)
  hist_axes_weight(myfit)
  hist_noise_weight(myfit)
  hist_free_scalex(myfit)
  hist_free_scaley(myfit)
  hist_axes_scalex(myfit)
  hist_axes_scaley(myfit)
  dev.off()
}

vis_traceplots <- function(myfit,name){
  pdf(paste(name, '_traceplots.pdf',sep=''))
  marg = 0.2;
  t1 = rstan::traceplot(myfit, pars='mix_weights[2]', inc_warmup = TRUE) + theme(legend.position='none',plot.margin=unit(c(marg,marg,marg,marg), "cm")) + ggtitle('Trace of mix weight') + xlab(NULL)+ylab(NULL)
  t2 = rstan::traceplot(myfit, pars='noise_weight', inc_warmup = TRUE) + theme(legend.position='none',plot.margin=unit(c(marg,marg,marg,marg), "cm")) + ggtitle('Trace of noise weight') + xlab(NULL)+ylab(NULL)
  t3 = rstan::traceplot(myfit, pars='xaligned_sd', inc_warmup = TRUE) + theme(legend.position='none',plot.margin=unit(c(marg,marg,marg,marg), "cm")) + ggtitle('Trace of axes x scale') + xlab(NULL)+ylab(NULL)
  t4 = rstan::traceplot(myfit, pars='yaligned_sd', inc_warmup = TRUE) + theme(legend.position='none',plot.margin=unit(c(marg,marg,marg,marg), "cm")) + ggtitle('Trace of axes y scale') + xlab(NULL)+ylab(NULL)
  t5 = rstan::traceplot(myfit, pars='axes_mix_weight', inc_warmup = TRUE) + theme(legend.position='none',plot.margin=unit(c(marg,marg,marg,marg), "cm")) + ggtitle('Trace of axes mix weight') + xlab(NULL)+ylab(NULL)
  t6 = rstan::traceplot(myfit, pars='free_scale[1]', inc_warmup = TRUE) + theme(legend.position='none',plot.margin=unit(c(marg,marg,marg,marg), "cm")) + ggtitle('Trace of free x scale') + xlab(NULL)+ylab(NULL)
  t7 = rstan::traceplot(myfit, pars='free_scale[2]', inc_warmup = TRUE) + theme(legend.position='none',plot.margin=unit(c(marg,marg,marg,marg), "cm")) + ggtitle('Trace of free y scale') + xlab(NULL)+ylab(NULL)
  t8 = rstan::traceplot(myfit, pars='free_cor[1,2]', inc_warmup = TRUE) + theme(legend.position='none',plot.margin=unit(c(marg,marg,marg,marg), "cm")) + ggtitle('Trace of correlation') + xlab(NULL)+ylab(NULL)
  grid.arrange(t1,t2,t3,t4,t5,t6,t7,t8,ncol=2)
  dev.off()
}
