library(rstan)
library(MASS) # For mvrnorm
library(plyr)
library(coda)
library(gridExtra)
library(ggplot2)
library(coefplot)
source("../Model/Visualizations.R")
dataDir = "../singleunitdata/";
figureDir = "../figures/";


## FIGURE 6

load(paste(dataDir, 'fiveKfits.RData', sep = ""))
for (i in 1:length(fits)) {
pdf(paste(figureDir, '6.',i,'.pdf',sep=''))
par(mfrow=c(2,2), mar=c(2,2.5,2,1))
plot_5cov_iters(fits[[i]][[8]],fits[[i]])
hist_mix_weight(fits[[i]][[8]])
hist_corr(fits[[i]][[8]])
hist_noise_weight(fits[[i]][[8]])

dev.off()
}

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

pdf(paste(figureDir, '7.pdf',sep=''), width=10)
multiplot(p1,p2,p3)
dev.off()

#Figure 8/SUPPLEMENT?
#a correlation vs naive correlation
#b correlation difference for correlation vs naive correlation on significant sets
d = as.data.frame(d)
numsets=length(fits) #adjust this line to make the numsets = the number of pairs of variables compaired (Tommy put in 20)
df = data.frame(matrix(nrow = numsets*2))
df$corr = c(d$median.correlation, d$naive.correlation)
df$low = c(d$low.correlation, d$low.naive.correlation)
df$high = c(d$high.correlation, d$high.naive.correlation)
df$name = c(paste(d$set.name), paste(d$set.name))
df$col = c(rep('Model mixed-tuning correlation', numsets), rep('Whole population correlation', numsets))

p = ggplot(data = df, aes(x = corr, y = name, color = col)) + 
  geom_point(position = position_dodgev(height = 0.5)) + # could not find function "position_dodgev"
  geom_errorbarh(aes(xmin = low, xmax = high), position = position_dodgev(height = 0.5))

ggsave(paste(figureDir, 'supp1.pdf', sep=''))

namesTable <- read.csv(paste(figureDir, 'namesTable.csv', sep = ""), header=T)
p = ggplot(as.data.frame(table(namesTable)), aes(x=gender, y = Freq, fill=fraud)) + geom_bar(stat="identity")


pdf(paste(figureDir, 'namesTable.pdf', sep=''), width=12)
grid.table(namesTable)
dev.off()
