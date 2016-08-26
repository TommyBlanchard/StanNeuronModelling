library(rstan)
library(MASS) # For mvrnorm
library(plyr)
library(coda)
library(gridExtra)
source("../Model/Visualizations.R")
dataDir = "../singleunitdata/";

#Figure 1 old methods
#a) beta x beta plot
#b) pie chart

make_data <- function(betas,N_CELLS,DATA_PER_CELL, noise){
  d <- NULL
  obs <- NULL
  pvals <- array(0,dim=c(N_CELLS,2))
  sigma <- array(0,dim=c(N_CELLS,2,2))
  for(r in 1:nrow(betas)) {
    b <- betas[r,]
    
    m <- matrix(runif(DATA_PER_CELL*2), ncol=2)
    
    err <- rnorm(DATA_PER_CELL, 0.0, noise)
    d <- data.frame(cell=r, x1=m[,1], x2=m[,2], y= m%*%b+err)
    fits <- lm(y ~ x1 + x2, data = d);
    obs <- rbind(obs,fits$coefficients[2:3]);
    fits.sum <- summary(fits)
    pvals[r,] <- coef(fits.sum)[2:3,4]; #pvalues for the coefficients
    sigma[r,,] <- vcov(fits)[2:3,2:3]; #covariance matrix for the coefficients
  }
  # Construct the data to send
  data <- list(N_CELLS=N_CELLS, 
               dim=2,
               sigma=sigma,
               alpha_cov_mix=array(rep(0.5,2),dim=2), # dirichlet prior. Need to declare as array to work with M=1
               X = obs,
               pvals = pvals
  )      
  return(data);
}
N_CELLS = 200;
DATA_PER_CELL = 50;
betas1 <- rbind(mvrnorm(n=N_CELLS, mu=c(0,0), Sigma=matrix(c(.1,0, 0,.1), nrow=2)))

data1= make_data(betas1,N_CELLS,DATA_PER_CELL, 0.4);
DATA_PER_CELL = 500;
data2= make_data(betas1,N_CELLS,DATA_PER_CELL, 0.4);

pdf(paste(dataDir, 'figures/Figure1.pdf',sep=''), width=8,height=8)
par(mfrow=c(2,2))
plot(data1$X, pch=20,ylim=c(-max(abs(c(data1$X))),max(abs(c(data1$X)))),xlim=c(-max(abs(c(data1$X))),max(abs(c(data1$X)))),cex.main=1,cex.axis=1.1, font.lab=2, font.main=2, font.axis=2, xlab = 'X1 tuning', ylab = 'X2 tuning')
abline(0,0,0,0)
ind = (data1$pvals[,1] >0.05) &  (data1$pvals[,2] >0.05);
mycol = col2rgb(palette()[2])/255
points(data1$X[ind,],pch=20, col=rgb(mycol[1],mycol[2],mycol[3]))

ind = xor(data1$pvals[,1] <0.05,data1$pvals[,2] <0.05);
mycol = col2rgb(palette()[3])/255
points(data1$X[ind,],pch=20, col=rgb(mycol[1],mycol[2],mycol[3]))

ind = (data1$pvals[,1] <0.05) &  (data1$pvals[,2] <0.05);
mycol = col2rgb(palette()[4])/255
points(data1$X[ind,],pch=20, col=rgb(mycol[1],mycol[2],mycol[3]))
title(main="A",adj=0)

labels = c('Neither', 'X1', 'X2', 'Both')
pie(c(sum((data1$pvals[,1] >0.05) &  (data1$pvals[,2] >0.05)),
      sum((data1$pvals[,1] <0.05) &  (data1$pvals[,2] >0.05)),
      sum((data1$pvals[,1] >0.05) &  (data1$pvals[,2] <0.05)),
      sum((data1$pvals[,1] <0.05) &  (data1$pvals[,2] <0.05))), 
    labels=labels, col=palette()[c(2:3,3,4)])
title(main="B",adj=0)

plot(data2$X, pch=20,ylim=c(-max(abs(c(data2$X))),max(abs(c(data2$X)))),xlim=c(-max(abs(c(data2$X))),max(abs(c(data2$X)))),cex.main=1,cex.axis=1.1, font.lab=2, font.main=2, font.axis=2, xlab = 'X1 tuning', ylab = 'X2 tuning')
abline(0,0,0,0)
ind = (data2$pvals[,1] >0.05) &  (data2$pvals[,2] >0.05);
mycol = col2rgb(palette()[2])/255
points(data2$X[ind,],pch=20, col=rgb(mycol[1],mycol[2],mycol[3]))

ind = xor(data2$pvals[,1] <0.05,data2$pvals[,2] <0.05);
mycol = col2rgb(palette()[3])/255
points(data2$X[ind,],pch=20, col=rgb(mycol[1],mycol[2],mycol[3]))

ind = (data2$pvals[,1] <0.05) &  (data2$pvals[,2] <0.05);
mycol = col2rgb(palette()[4])/255
points(data2$X[ind,],pch=20, col=rgb(mycol[1],mycol[2],mycol[3]))
title(main="C",adj=0)

labels = c('Neither', 'X1', 'X2', 'Both')
pie(c(sum((data2$pvals[,1] >0.05) &  (data2$pvals[,2] >0.05)),
      sum((data2$pvals[,1] <0.05) &  (data2$pvals[,2] >0.05)),
      sum((data2$pvals[,1] >0.05) &  (data2$pvals[,2] <0.05)),
      sum((data2$pvals[,1] <0.05) &  (data2$pvals[,2] <0.05))), 
    labels=labels, col=palette()[c(2:3,3,4)])
title(main="D",adj=0)
dev.off()

#Figure 2 model & simulated data
#a) model diagram

#b) mixed selective dists
modelcode <- paste(readLines('model.stan'), collapse = '\n')
ITER = 1000
CHAINS = 4
fit1 <- stan(model_code=modelcode, data=data1, iter=ITER, chains=CHAINS);

pdf(paste(dataDir, 'figures/2b.pdf',sep=''), width=4,height=7.5)
layout(matrix(c(1,2,3,4,5,5,5,5), 4, byrow='true'), widths = rep.int(1, 2),
       heights = rep.int(1, 4), respect = TRUE)
par(pty="s",
    mar = c(2, 1, 1, 0) + 0.1,
    oma = c(2,2,0,0) + 0.1)
plot_data_covs(data1)
plot_5free_iters(fit1,data1)
plot_5axes_iters(fit1,data1)
plot_5noise_iters(fit1,data1)
plot_5cov_iters(fit1,data1)
title(xlab = 'X1 tuning',
      ylab = 'X2 tuning',
      outer = TRUE,
      line=0.1,
      font.lab=2,
      cex.lab = 1.5)

dev.off()


#c)

betas2 <- rbind(mvrnorm(n=(N_CELLS/4), mu=c(0,0), Sigma=matrix(c(0,0, 0,0), nrow=2)),
                mvrnorm(n=(N_CELLS/4), mu=c(0,0),Sigma=matrix(c(.1,0, 0,0), nrow=2)),
                mvrnorm(n=(N_CELLS/4), mu=c(0,0), Sigma=matrix(c(0,0, 0,.1), nrow=2)),
                mvrnorm(n=(N_CELLS/4), mu=c(0,0), Sigma=matrix(c(.1,0, 0,.1), nrow=2)))     
data2= make_data(betas2,N_CELLS,DATA_PER_CELL, 0.4);
fit2 <- stan(model_code=modelcode, data=data2, iter=ITER, chains=CHAINS);

pdf(paste(dataDir, 'figures/2c.pdf',sep=''), width=4,height=7.5)
layout(matrix(c(1,2,3,4,5,5,5,5), 4, byrow='true'), widths = rep.int(1, 2),
       heights = rep.int(1, 4), respect = TRUE)
par(pty="s",
    mar = c(2, 1, 1, 0) + 0.1,
    oma = c(2,2,0,0) + 0.1)
plot_data_covs(data2)
plot_5free_iters(fit2,data2)
plot_5axes_iters(fit2,data2)
plot_5noise_iters(fit2,data2)
plot_5cov_iters(fit2,data2)
title(xlab = 'X1 tuning',
      ylab = 'X2 tuning',
      outer = TRUE,
      line=0.1,
      font.lab=2,
      cex.lab = 1.5)
dev.off()

#Figure 3 posteriors for models in fig2

pdf(paste(dataDir, 'figures/3.pdf',sep=''))
par(mfrow=c(2,3), mar=c(2,2.5,2,1))
hist_mix_weight(fit1)
hist_axes_weight(fit1)
hist_noise_weight(fit1)
hist_mix_weight(fit2)
hist_axes_weight(fit2)
hist_noise_weight(fit2)
dev.off()

#Figure 4 Correlated tuning
#d&e) correlated mixed selectivity
betas3 <- rbind(mvrnorm(n=N_CELLS, mu=c(0,0), Sigma=matrix(c(.1,.05, .05,.1), nrow=2)))
data3= make_data(betas3,N_CELLS,DATA_PER_CELL, 0.4);
fit3 <- stan(model_code=modelcode, data=data3, iter=ITER, chains=CHAINS);
betas4 <- rbind(mvrnorm(n=(N_CELLS/4), mu=c(0,0), Sigma=matrix(c(0,0, 0,0), nrow=2)),
                mvrnorm(n=(N_CELLS/4), mu=c(0,0),Sigma=matrix(c(.1,0, 0,0), nrow=2)),
                mvrnorm(n=(N_CELLS/4), mu=c(0,0), Sigma=matrix(c(0,0, 0,.1), nrow=2)),
                mvrnorm(n=(N_CELLS/4), mu=c(0,0), Sigma=matrix(c(.1,.05, .05,.1), nrow=2)))     
data4= make_data(betas4,N_CELLS,DATA_PER_CELL, 0.4);
fit4 <- stan(model_code=modelcode, data=data4, iter=ITER, chains=CHAINS);

pdf(paste(dataDir, 'figures/4a.pdf',sep=''), width=4,height=4)
plot_5cov_iters(fit3,data3)
dev.off()

pdf(paste(dataDir, 'figures/4b.pdf',sep=''), width=4,height=4)
plot_5cov_iters(fit4,data4)
dev.off()

pdf(paste(dataDir, 'figures/4.pdf',sep=''))
par(mfrow=c(4,2), mar=c(2,2.5,2,1))
hist_mix_weight(fit3)
hist_corr(fit3)
hist_axes_weight(fit3)
hist_noise_weight(fit3)
hist_mix_weight(fit4)
hist_corr(fit4)
hist_axes_weight(fit4)
hist_noise_weight(fit4)
dev.off()

#Figure 5: behavior when data limited
#a) mixed-selective with noise
#b) pure-tuned with noise
#c) mix of all three

## FIGURE 6

load(paste(dataDir, 'fiveKfits.RData', sep = ""))

pdf(paste(dataDir, 'figures/6.pdf',sep=''))
par(mfrow=c(2,2), mar=c(2,2.5,2,1))
plot_5cov_iters(fits$ofcCuriosity_Rew_Info[[8]],fits$ofcCuriosity_Rew_Info)
hist_mix_weight(fits$ofcCuriosity_Rew_Info[[8]])
hist_corr(fits$ofcCuriosity_Rew_Info[[8]])
hist_noise_weight(fits$ofcCuriosity_Rew_Info[[8]])

dev.off()

## FIGURE 7
#a) mixed-weight proporiton of signal weight
#b) pure-weight as proportion of total weight
d <- read.csv(paste(dataDir, 'weightsTable.csv', sep = ""), header=T)
p1 <- ggplot(d, aes(y=X, x=median.mixed.tuning.signal.weight, xmin=low.mixed.tuning.signal.weight, xmax=high.mixed.tuning.signal.weight)) + 
  geom_point() + 
  geom_errorbarh() +
  theme_bw() +
  xlim(0,1) +
  labs(x = "Mixed-tuning weight as proportion of signal weight", y='Data set')+
  scale_y_continuous(breaks = 1:10*2)

p2 <- ggplot(d, aes(y=X, x=median.pure.weight, xmin=low.pure.weight, xmax=high.pure.weight)) + 
  geom_point() + 
  geom_errorbarh() +
  theme_bw() +
  xlim(0,1) +
  labs(x = "Pure-tuning weight as proportion of total weight", y='Data set')+
  scale_y_continuous(breaks = 1:10*2)

p3 <- ggplot(d, aes(y=X, x=median.noise.weight, xmin=low.noise.weight, xmax=high.noise.weight)) + 
  geom_point() + 
  geom_errorbarh() +
  theme_bw() +
  xlim(0,1) +
  labs(x = "No-tuning weight as proportion of total weight", y='Data set')+
  scale_y_continuous(breaks = 1:10*2)

pdf(paste(dataDir, 'figures/7.pdf',sep=''), width=10)
multiplot(p1,p2,p3)
dev.off()

#Figure 8/SUPPLEMENT?
#a correlation vs naive correlation
#b correlation difference for correlation vs naive correlation on significant sets
df = data.frame(matrix(nrow = 40))
d = as.data.frame(d)
df$corr = c(d$median.correlation, d$naive.correlation)
df$low = c(d$low.correlation, d$low.naive.correlation)
df$high = c(d$high.correlation, d$high.naive.correlation)
df$name = c(paste(d$set.name), paste(d$set.name))
df$col = c(rep('Model mixed-tuning correlation', 20), rep('Whole population correlation', 20))

p = ggplot(data = df, aes(x = corr, y = name, color = col)) + 
  geom_point(position = position_dodgev(height = 0.5)) +
  geom_errorbarh(aes(xmin = low, xmax = high), position = position_dodgev(height = 0.5))

ggsave(paste(dataDir, 'figures/supp1.pdf', sep=''))

namesTable <- read.csv(paste(dataDir, 'figures/namesTable.csv', sep = ""), header=T)
p = ggplot(as.data.frame(table(namesTable)), aes(x=gender, y = Freq, fill=fraud)) + geom_bar(stat="identity")


pdf(paste(dataDir, 'figures/namesTable.pdf', sep=''), width=12)
grid.table(namesTable)
dev.off()
