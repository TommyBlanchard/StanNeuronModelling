makeTableRow <-function(fit){
  pars = extract(fit[[8]])
  free_weight = pars$mix_weights[,2]
  free_weight_high= HPDinterval(as.mcmc(free_weight), prob=0.95)[2];
  free_weight_low = HPDinterval(as.mcmc(free_weight), prob=0.95)[1];
  
  noise_weight = array(pars$noise_weight)
  noise_weight_high= HPDinterval(as.mcmc(noise_weight), prob=0.95)[2];
  noise_weight_low = HPDinterval(as.mcmc(noise_weight), prob=0.95)[1];
  
  total_pure_weight = (1 - free_weight)*(1 - noise_weight)
  pure_weight_high= HPDinterval(as.mcmc(total_pure_weight), prob=0.95)[2];
  pure_weight_low = HPDinterval(as.mcmc(total_pure_weight), prob=0.95)[1];
  
  free_cor = pars$free_cor[,1,2]
  free_cor_high = HPDinterval(as.mcmc(free_cor), prob=0.95)[2];
  free_cor_low = HPDinterval(as.mcmc(free_cor), prob=0.95)[1];
  
  naive_cor = cor.test(fit$X[,1],fit$X[,2])$estimate
  naive_cor_low = cor.test(fit$X[,1],fit$X[,2])$conf.int[1]
  naive_cor_high = cor.test(fit$X[,1],fit$X[,2])$conf.int[2]
  
  c(median(free_weight),free_weight_low,free_weight_high,median(total_pure_weight),pure_weight_low,pure_weight_high,median(noise_weight),noise_weight_low,noise_weight_high,median(free_cor),free_cor_low,free_cor_high, naive_cor, naive_cor_low, naive_cor_high)
}

table = NULL
for (i in 1:length(fits)) {
  table <- rbind(table,makeTableRow(fits[[i]]));
}
table <- cbind(levels(regdata$setName),table)
table[,1] = levels(regdata$setName)[c(1:6, 8, 10:14, 16:23)] #Took out 3 datasets that weren't converging as quick as the rest, temporary hack to get rid of those
colnames(table) <- c("set name", "median mixed-tuning signal weight", "low mixed-tuning signal weight", "high mixed-tuning signal weight", "median pure weight", "low pure weight", "high pure weight", "median noise weight", "low noise weight", "high noise weight",  "median correlation", "low correlation", "high correlation", "naive correlation", "low naive correlation", "high naive correlation")

write.csv(table,paste(dataDir, 'weightsTable.csv', sep = ""))

#check rhats
Rhat <- function(fit) {
  summary(fit)$summary[,"Rhat"]
}
rhats = NULL
for (i in 1:length(fits)) {
  rhats <- rbind(rhats,Rhat(fits[[i]][[8]]));
}

