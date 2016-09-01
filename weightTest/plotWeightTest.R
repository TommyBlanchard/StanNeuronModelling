library(ggplot2)

pSignal = NULL
pNoise = NULL
for (j in 1:3){
  dat = NULL
  for (i in 1:20){
    load(paste("weights/weights", i, ".", j, ".RData", sep=""))
    dat = rbind(dat, d)
  }
  if(j == 1)
    main='Mixed-Tuned Population'
  else if(j == 2)
    main='Pure-Tuned Population'
  else if(j == 3)
    main='Pure and Mixed Population'
  pSignal[[j]] <- ggplot(as.data.frame(dat), aes(y=free_weight_med, x=neurons, ymin=free_weight_low, ymax=free_weight_high)) + 
    geom_point() + 
    geom_errorbar() +
    theme_bw() +
    ylim(0,1) +
    labs(y = "Mixed-tuning weight", title=main)
  pNoise[[j]] <- ggplot(as.data.frame(dat), aes(y=noise_weight_med, x=neurons, ymin=noise_weight_low, ymax=noise_weight_high)) + 
    geom_point() + 
    geom_errorbar() +
    theme_bw() +
    ylim(0,1) +
    labs(y = "No-tuning weight", title=main)
}

pdf(paste(dataDir, 'figures/5.pdf',sep=''), width=8,height=8)
multiplot(plotlist=c(pSignal,pNoise),cols=2)
dev.off()
