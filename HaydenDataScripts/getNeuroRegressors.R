library(plyr)

dataDir = "/StanNeuronModelling/singleunitdata/";

runReg = function(x) {
  tryCatch({
  if(any(x$regInd ==1)){
    print(x$setName[1])
    print(x$neuron[1])
    #If regInd has 2 paste, we do 2 different regressions. 
    #(e.g. spikes comes from different epochs/trials for the two variables being tested)
    fits1 <- lm((scale(as.double(paste(spikes)))) ~ scale(as.double(paste(X1))), data = x[which(x$regInd == 0),])
    fits2 <- lm(scale(as.double(paste(spikes))) ~ scale(as.double(paste(X2))), data = x[which(x$regInd == 1),])
    beta <- c(fits1$coefficients[2],fits2$coefficients[2]);
    sigma <- c(vcov(fits1)[2,2], 0, 0, vcov(fits2)[2,2]) #no covariance if spike counts are taken from separate periods
  }else{
    #If 'spikes' comes from the same epoch for both variables, regInd should only have 1 level and everything should be run as one regression
    fits <- lm(scale(as.double(paste(spikes))) ~ scale(as.double(paste(X1))) + scale(as.double(paste(X2))), data = x);
    beta <- fits$coefficients[2:3];
    sigma <- vcov(fits)[2:3,2:3]; #covariance matrix for the coefficients
  }
  dim(beta) <- NULL
  dim(sigma) <- NULL
  
  setNames(c(beta, sigma), c("beta1", "beta2", "sigma11", "sigma21", "sigma12", "sigma22"))
  },error = function(e){
    #lm throws an error if there wasn't a single spike in any trial for a neuron, this catches that error and just returns na
    print(paste('error with neuron ',x$neuron[1]))
  setNames(c(NA,NA,NA,NA,NA,NA), c("beta1", "beta2", "sigma11", "sigma21", "sigma12", "sigma22"))
  })
}

######### OFC CURIOSITY
ofcCuriosity = read.csv(paste(dataDir, "ofcCuriosity.csv", sep = ""))
nrows = dim(ofcCuriosity)[1]
setName = rep('ofcCuriosity_Rew_Info', nrows)

#Rew vs Info, (FIRST EPOCH)
firstRew = ofcCuriosity$rightrew
firstRew[ofcCuriosity$leftfirst == 1] = ofcCuriosity$leftrew[ofcCuriosity$leftfirst == 1]
firstType = ofcCuriosity$righttype
firstType[ofcCuriosity$leftfirst == 1] = ofcCuriosity$lefttype[ofcCuriosity$leftfirst == 1]

X1 = firstRew
X2 = firstType

spikes = rowSums(ofcCuriosity[,grep("^offerPsth263$", colnames(ofcCuriosity)):grep("^offerPsth287$", colnames(ofcCuriosity))])
neuron = ofcCuriosity$neuron
#regInd is 0 if all data is going in the same regression, otherwise indicated 0 or 1 for which of two regressions
#this data point should be used for
regInd = rep(0,nrows) 
data = cbind(setName, neuron, spikes, X1, X2, as.integer(regInd))

data = as.data.frame(data)

#RewInfo vs RewNoInfo, (FIRST EPOCH)
setName = rep('ofcCuriosity_RewInfo_RewNoInfo', nrows)
X1= firstRew
X1[firstType == 2] = 0
X2= firstRew
X2[firstType == 1] = 0 
spikes = rowSums(ofcCuriosity[,grep("^offerPsth263$", colnames(ofcCuriosity)):grep("^offerPsth287$", colnames(ofcCuriosity))])

neuron = ofcCuriosity$neuron
regInd = rep(0,nrows)

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#Opt1 rew vs Opt2 rew (SECOND EPOCH)
setName = rep('ofcCuriosity_opt1Rew_opt2Rew', nrows)
secondRew = ofcCuriosity$rightrew
secondRew[ofcCuriosity$leftfirst == 0] = ofcCuriosity$leftrew[ofcCuriosity$leftfirst == 0]
X1= firstRew
X2= secondRew
spikes = rowSums(ofcCuriosity[,grep("^offerPsth301$", colnames(ofcCuriosity)):grep("^offerPsth325$", colnames(ofcCuriosity))])

neuron = ofcCuriosity$neuron
regInd = rep(0,nrows)

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#left rew vs right rew (FIRST & SECOND EPOCH)
setName = rep('ofcCuriosity_leftRew_rightRew', nrows)
X1= rep(ofcCuriosity$rightrew, 2)
X2= rep(ofcCuriosity$leftrew, 2)

regInd = c(ofcCuriosity$leftFirst, !ofcCuriosity$leftfirst)

spikes = c(rowSums(ofcCuriosity[,grep("^offerPsth263$", colnames(ofcCuriosity)):grep("^offerPsth287$", colnames(ofcCuriosity))]),
           rowSums(ofcCuriosity[,grep("^offerPsth301$", colnames(ofcCuriosity)):grep("^offerPsth325$", colnames(ofcCuriosity))]))

neuron = rep(ofcCuriosity$neuron,2)

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#first offer rew vs chosen rew (TWO EPOCHS/REGRESSIONS - first epoch and choice epoch)
setName = rep('ofcCuriosity_firstRew_chosenRew', nrows*2)
chosenRew = ofcCuriosity$rightrew
chosenRew[ofcCuriosity$leftaccept == 1] = ofcCuriosity$leftrew[ofcCuriosity$leftaccept == 1]
X1= c(firstRew,rep(0,nrows))
X2= c(rep(0,nrows),chosenRew)
spikes = c(rowSums(ofcCuriosity[,grep("^offerPsth263$", colnames(ofcCuriosity)):grep("^offerPsth287$", colnames(ofcCuriosity))]),
           rowSums(ofcCuriosity[,grep("^offerPsth338$", colnames(ofcCuriosity)):grep("^offerPsth362$", colnames(ofcCuriosity))]))

neuron = rep(ofcCuriosity$neuron,2)
regInd = c(rep(0,nrows), rep(1,nrows))

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#first offer rew vs outcome rew (TWO EPOCHS/REGRESSIONS - first epoch and outcome epoch)
setName = rep('ofcCuriosity_firstRew_outcomeRew', nrows*2)
X1= c(firstRew,rep(0,nrows))
X2= c(rep(0,nrows),ofcCuriosity$reward)
spikes = c(rowSums(ofcCuriosity[,grep("^offerPsth263$", colnames(ofcCuriosity)):grep("^offerPsth287$", colnames(ofcCuriosity))]),
           rowSums(ofcCuriosity[,grep("^outcomePsth251$", colnames(ofcCuriosity)):grep("^offerPsth275$", colnames(ofcCuriosity))]))

neuron = rep(ofcCuriosity$neuron,2)
regInd = c(rep(0,nrows), rep(1,nrows))

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

######### OFC SIMULATION
ofcSimulation = read.csv(paste(dataDir, "ofcSimulation.csv", sep = ""))
neuron = ofcSimulation$neuron
nrows = dim(ofcSimulation)[1]

offerRew = c(75,100,150,200,250)[ofcSimulation$offerNum]
standardRew = c(150,175,200)[ofcSimulation$standardNum/10]

#experienced vs described offers
#setName = rep('ofcSimulation_exp_desc', nrows)
#X1 = offerRew
#X2 = offerRew
#spikes = rowSums(ofcSimulation[,grep("^offerPsth261$", colnames(ofcSimulation)):grep("^offerPsth275$", colnames(ofcSimulation))])
#regInd = ofcSimulation$experienceTrial
#data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#experienced offer rew vs outcome rew
choseStandard = (ofcSimulation$leftStandard & ofcSimulation$choseLeft) | (!ofcSimulation$leftStandard & !ofcSimulation$choseLeft)
chosenRew = offerRew
chosenRew[choseStandard] = standardRew[choseStandard]

setName = rep('ofcSimulation_exp_outcome', sum(ofcSimulation$experienceTrial & !choseStandard) + sum(ofcSimulation$experienceTrial & !choseStandard))

X1 = c(chosenRew[!choseStandard  & ofcSimulation$experienceTrial],rep(0,sum(!choseStandard  & ofcSimulation$experienceTrial)))
X2 = c(rep(0,sum(!choseStandard & ofcSimulation$experienceTrial == 1) ),offerRew[!choseStandard  & ofcSimulation$experienceTrial])
spikes = c(rowSums(ofcSimulation[!choseStandard & ofcSimulation$experienceTrial == 1,grep("^outcomePsth251$", colnames(ofcSimulation)):grep("^outcomePsth288$", colnames(ofcSimulation))]),
           rowSums(ofcSimulation[!choseStandard & ofcSimulation$experienceTrial == 1,grep("^offerPsth261$", colnames(ofcSimulation)):grep("^offerPsth288$", colnames(ofcSimulation))]))

neuron = c(ofcSimulation$neuron[ofcSimulation$experienceTrial == 1 & !choseStandard],ofcSimulation$neuron[!choseStandard  & ofcSimulation$experienceTrial])
regInd = c(rep(0,sum(!choseStandard & ofcSimulation$experienceTrial == 1)),rep(1,sum(!choseStandard  & ofcSimulation$experienceTrial)))
data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#described offer rew vs outcome rew
#setName = rep('ofcSimulation_desc_outcome', sum(!ofcSimulation$experienceTrial & !choseStandard) + sum(!ofcSimulation$experienceTrial & !choseStandard))

#X1 = c(chosenRew[!choseStandard  & !ofcSimulation$experienceTrial],rep(0,sum(!choseStandard  & !ofcSimulation$experienceTrial)))
#X2 = c(rep(0,sum(!choseStandard & !ofcSimulation$experienceTrial == 1) ),offerRew[!choseStandard  & !ofcSimulation$experienceTrial])
#spikes = c(rowSums(ofcSimulation[!choseStandard & !ofcSimulation$experienceTrial == 1,grep("^outcomePsth251$", colnames(ofcSimulation)):grep("^outcomePsth288$", colnames(ofcSimulation))]),
#           rowSums(ofcSimulation[!choseStandard & !ofcSimulation$experienceTrial == 1,grep("^offerPsth261$", colnames(ofcSimulation)):grep("^offerPsth288$", colnames(ofcSimulation))]))

#neuron = c(ofcSimulation$neuron[!ofcSimulation$experienceTrial == 1 & !choseStandard],ofcSimulation$neuron[!choseStandard  & !ofcSimulation$experienceTrial])
#regInd = c(rep(0,sum(!choseStandard & !ofcSimulation$experienceTrial == 1)),rep(1,sum(!choseStandard  & !ofcSimulation$experienceTrial)))
#data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

######### dACC, DIET SELECTION 
daccDietSelection = read.csv(paste(dataDir, "daccDietSelection.csv", sep = ""))
neuron = daccDietSelection$neuron
nrows = dim(daccDietSelection)[1]

#rew vs delay, (500ms after option appearance)
setName = rep('daccDietSelection_Rew_Delay', nrows)
X1= daccDietSelection$offerRewSize
X2= daccDietSelection$offerDelay
spikes = rowSums(daccDietSelection[,grep("^offerPsth276$", colnames(daccDietSelection)):grep("^offerPsth300$", colnames(daccDietSelection))])
regInd = rep(0,nrows)
data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#offer rew vs outcome rew (stim appearance, outcome epoch)
setName = rep('daccDietSelection_Offer_Outcome', nrows*2)
X1= c(daccDietSelection$offerRewSize,rep(0,nrows))
X2= c(rep(0,nrows),daccDietSelection$offerRewSize)
spikes = c(rowSums(daccDietSelection[,grep("^offerPsth276$", colnames(daccDietSelection)):grep("^offerPsth300$", colnames(daccDietSelection))]),
           rowSums(daccDietSelection[,grep("^outcomePsth251$", colnames(daccDietSelection)):grep("^offerPsth275$", colnames(daccDietSelection))]))

neuron = rep(daccDietSelection$neuron,2)
regInd = c(rep(0,nrows), rep(1,nrows))

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

######### dACC, TOKENS
daccTokens = read.csv(paste(dataDir, "daccTokens.csv", sep = ""))
#Note: the data is in spikes/second instead of spikes/bin, shouldn't matter since we're zscoring
#Just get rid of trials where an option's top value == bottom value (equivalent to a safe option)
daccTokens <- daccTokens[daccTokens$leftBottomTokenValue != daccTokens$leftTopTokenValue & daccTokens$rightBottomTokenValue != daccTokens$rightTopTokenValue,]

#Also get rid of the trials with nan values for top prob
daccTokens <- daccTokens[!is.na(daccTokens$rightTopProb),]
neuron = daccTokens$neuron
nrows = dim(daccTokens)[1]

#stakes vs probability (epoch 1)
setName = rep('daccTokens_Stakes_Probability', nrows)

leftMax = apply(daccTokens[,c(grep("^leftTopTokenValue$", colnames(daccTokens)),grep("^leftBottomTokenValue$", colnames(daccTokens)))],1,max)
rightMax = apply(daccTokens[,c(grep("^rightTopTokenValue$", colnames(daccTokens)),grep("^rightBottomTokenValue$", colnames(daccTokens)))],1,max)
leftMaxInd = apply(daccTokens[,c(grep("^leftTopTokenValue$", colnames(daccTokens)),grep("^leftBottomTokenValue$", colnames(daccTokens)))],1,which.max)
rightMaxInd = apply(daccTokens[,c(grep("^rightTopTokenValue$", colnames(daccTokens)),grep("^rightBottomTokenValue$", colnames(daccTokens)))],1,which.max)
leftMaxProb = daccTokens$leftTopProb
leftMaxProb[as.integer(leftMaxInd) == 2] = 100 - leftMaxProb[as.integer(leftMaxInd) == 2]
rightMaxProb = daccTokens$rightTopProb
rightMaxProb[as.integer(rightMaxInd) == 2] = 100 - daccTokens$rightTopProb[as.integer(rightMaxInd) == 2]

firstMax = leftMax
firstMax[daccTokens$orderOfLeft == 2] = rightMax[daccTokens$orderOfLeft == 2]
firstMaxProb = leftMaxProb
firstMaxProb[daccTokens$orderOfLeft == 2] = rightMaxProb[daccTokens$orderOfLeft == 2]

X1 = firstMax
X2 = firstMaxProb
spikes = rowSums(daccTokens[,grep("^offerPsth251$", colnames(daccTokens)):grep("^offerPsth275$", colnames(daccTokens))])
regInd = rep(0,nrows)
data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#expected value option 1, value option 2
setName = rep('daccTokens_ev1_ev2', nrows)

leftEV = daccTokens$leftTopTokenValue*daccTokens$leftTopProb/100 + daccTokens$leftBottomTokenValue*(100 - daccTokens$leftTopProb)/100
rightEV = daccTokens$rightTopTokenValue*daccTokens$rightTopProb/100 + daccTokens$rightBottomTokenValue*(100 - daccTokens$rightTopProb)/100
firstEV = leftEV
firstEV[daccTokens$orderOfLeft == 2] = rightEV[daccTokens$orderOfLeft == 2]
secondEV = leftEV
secondEV[daccTokens$orderOfLeft == 1] = rightEV[daccTokens$orderOfLeft == 1]

X1 = firstEV
X2 = secondEV
spikes = rowSums(daccTokens[,grep("^offerPsth288$", colnames(daccTokens)):grep("^offerPsth312$", colnames(daccTokens))])
regInd = rep(0,nrows)
data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#left offer expected value vs right offer expected value (first and second epochs)
setName = rep('daccTokens_leftEV_rightEV', nrows)

X1 = rep(leftEV,2)
X2 = rep(rightEV,2)

spikes = c(rowSums(daccTokens[,grep("^offerPsth251$", colnames(daccTokens)):grep("^offerPsth275$", colnames(daccTokens))]),
           rowSums(daccTokens[,grep("^offerPsth288$", colnames(daccTokens)):grep("^offerPsth312$", colnames(daccTokens))]))

regInd = c(daccTokens$orderOfLeft == 2, daccTokens$orderOfLeft == 1)

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#chosen value vs offered value (epoch 1 and choice epoch)
#setName = rep('daccTokens_firstEV_chosenEV', nrows*2)

#chosenEV = firstEV
#chosenEV[daccTokens$chosenOpt == 2] = secondEV[daccTokens$chosenOpt == 2]

#X1= c(firstEV,rep(0,nrows))
#X2= c(rep(0,nrows),chosenEV)
#spikes = c(rowSums(daccTokens[,grep("^offerPsth251$", colnames(daccTokens)):grep("^offerPsth275$", colnames(daccTokens))]),
#           rowSums(daccTokens[,grep("^offerPsth331$", colnames(daccTokens)):grep("^offerPsth355$", colnames(daccTokens))]))

#neuron = rep(daccTokens$neuron,2)
#regInd = c(rep(0,nrows), rep(1,nrows))

#data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

######### VS, gambling:
vsGambling = read.csv(paste(dataDir, "VSGambling.csv", sep = ""))

#remove invalid & safe trials
vsGambling <- vsGambling[vsGambling$validTrial==1,]
vsGambling <- vsGambling[(vsGambling$opt1Prob < 1),]
vsGambling <- vsGambling[(vsGambling$opt2Prob < 1),]
neuron = vsGambling$neuron
nrows = dim(vsGambling)[1]

#reward vs probability (epoch 1)
setName = rep('vsGambling_rew_prob', nrows)
X1= vsGambling$opt1Rew
X2= vsGambling$opt1Prob
spikes = rowSums(vsGambling[,grep("^offerPsth251$", colnames(vsGambling)):grep("^offerPsth275$", colnames(vsGambling))])
regInd = rep(0,nrows)
data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#value option 1, value option 2 (epoch 2)
setName = rep('vsGambling_ev1_ev2', nrows)
X1= vsGambling$opt1EV
X2= vsGambling$opt2EV
spikes = rowSums(vsGambling[,grep("^offerPsth301$", colnames(vsGambling)):grep("^offerPsth325$", colnames(vsGambling))])
regInd = rep(0,nrows)
data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#left rew vs right rew (FIRST & SECOND EPOCH)
setName = rep('vsGambling_leftEV_rightEV', nrows)
leftEV = vsGambling$opt1EV
leftEV[vsGambling$opt1Side != 1] = vsGambling$opt2EV[vsGambling$opt1Side != 1]
rightEV = vsGambling$opt2EV
rightEV[vsGambling$opt1Side != 1] = vsGambling$opt1EV[vsGambling$opt1Side != 1]

X1= rep(leftEV, 2)
X2= rep(rightEV, 2)

regInd = c(vsGambling$opt1Side == 2, vsGambling$opt1Side == 1)

spikes = c(rowSums(vsGambling[,grep("^offerPsth251$", colnames(vsGambling)):grep("^offerPsth275$", colnames(vsGambling))]),
           rowSums(vsGambling[,grep("^offerPsth301$", colnames(vsGambling)):grep("^offerPsth325$", colnames(vsGambling))]))

neuron = rep(vsGambling$neuron,2)

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#chosen value vs offered value (epoch 1 and choice epoch)
setName = rep('vsGambling_firstEV_chosenEV', nrows*2)
chosenEV = vsGambling$opt1EV
chosenEV[vsGambling$chosenOpt == 2] = vsGambling$opt2EV[vsGambling$chosenOpt == 2]

X1= c(vsGambling$opt1EV,rep(0,nrows))
X2= c(rep(0,nrows),chosenEV)
spikes = c(rowSums(vsGambling[,grep("^offerPsth251$", colnames(vsGambling)):grep("^offerPsth275$", colnames(vsGambling))]),
           rowSums(vsGambling[,grep("^offerPsth331$", colnames(vsGambling)):grep("^offerPsth356$", colnames(vsGambling))]))

neuron = rep(vsGambling$neuron,2)
regInd = c(rep(0,nrows), rep(1,nrows))

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

######### vmpfc, gambling:
vmpfcGambling = read.csv(paste(dataDir, "vmpfcGambling.csv", sep = ""))

#remove invalid & safe trials
vmpfcGambling <- vmpfcGambling[vmpfcGambling$validTrial==1,]
vmpfcGambling <- vmpfcGambling[(vmpfcGambling$opt1Prob < 1),]
vmpfcGambling <- vmpfcGambling[(vmpfcGambling$opt2Prob < 1),]
neuron = vmpfcGambling$neuron
nrows = dim(vmpfcGambling)[1]

#reward vs probability (epoch 1)
setName = rep('vmpfcGambling_rew_prob', nrows)
X1= vmpfcGambling$opt1Rew
X2= vmpfcGambling$opt1Prob
spikes = rowSums(vmpfcGambling[,grep("^offerPsth251$", colnames(vmpfcGambling)):grep("^offerPsth275$", colnames(vmpfcGambling))])
regInd = rep(0,nrows)
data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#value option 1, value option 2 (epoch 2)
setName = rep('vmpfcGambling_ev1_ev2', nrows)
X1= vmpfcGambling$opt1EV
X2= vmpfcGambling$opt2EV
spikes = rowSums(vmpfcGambling[,grep("^offerPsth301$", colnames(vmpfcGambling)):grep("^offerPsth325$", colnames(vmpfcGambling))])
regInd = rep(0,nrows)
data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#left rew vs right rew (FIRST & SECOND EPOCH)
setName = rep('vmpfcGambling_leftEV_rightEV', nrows)
leftEV = vmpfcGambling$opt1EV
leftEV[vmpfcGambling$opt1Side != 1] = vmpfcGambling$opt2EV[vmpfcGambling$opt1Side != 1]
rightEV = vmpfcGambling$opt2EV
rightEV[vmpfcGambling$opt1Side != 1] = vmpfcGambling$opt1EV[vmpfcGambling$opt1Side != 1]

X1= rep(leftEV, 2)
X2= rep(rightEV, 2)

regInd = c(vmpfcGambling$opt1Side == 2, vmpfcGambling$opt1Side == 1)

spikes = c(rowSums(vmpfcGambling[,grep("^offerPsth251$", colnames(vmpfcGambling)):grep("^offerPsth275$", colnames(vmpfcGambling))]),
           rowSums(vmpfcGambling[,grep("^offerPsth301$", colnames(vmpfcGambling)):grep("^offerPsth325$", colnames(vmpfcGambling))]))

neuron = rep(vmpfcGambling$neuron,2)

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

#chosen value vs offered value (epoch 1 and choice epoch)
setName = rep('vmpfcGambling_firstEV_chosenEV', nrows*2)
chosenEV = vmpfcGambling$opt1EV
chosenEV[vmpfcGambling$chosenOpt == 2] = vmpfcGambling$opt2EV[vmpfcGambling$chosenOpt == 2]

X1= c(vmpfcGambling$opt1EV,rep(0,nrows))
X2= c(rep(0,nrows),chosenEV)
spikes = c(rowSums(vmpfcGambling[,grep("^offerPsth251$", colnames(vmpfcGambling)):grep("^offerPsth275$", colnames(vmpfcGambling))]),
           rowSums(vmpfcGambling[,grep("^offerPsth331$", colnames(vmpfcGambling)):grep("^offerPsth356$", colnames(vmpfcGambling))]))

neuron = rep(vmpfcGambling$neuron,2)
regInd = c(rep(0,nrows), rep(1,nrows))

data = rbind(data, cbind(setName, neuron, spikes, X1, X2, as.integer(regInd)))

save(data,file=paste(dataDir, 'regressorInput.RData', sep = ""))

## RUN THE REGRESSION ON EVERYTHING  
names(data)[6] <- 'regInd'
#It's unclear why this needs to be done, but otherwise ddply throws an error
names(data$setName) <-NULL
names(data$neuron) <-NULL
names(data$spikes) <-NULL
names(data$X1) <-NULL
names(data$X2) <-NULL
names(data$regInd) <-NULL
regdata = ddply(data, c("setName", "neuron"), runReg)

#remove any nans caused by bad cells
regdata <- regdata[complete.cases(regdata),]

save(regdata,file=paste(dataDir, 'regdata.RData', sep = ""))


## SOME DEBUGGING TOOLS
# obs <- NULL
# ps <- NULL
# for (i in 1:121) {
#   fits <- lm(scale(as.double(paste(spikes))) ~ scale(as.double(paste(X1))) + scale(as.double(paste(X2))), data = data1[which(data1$neuron==i & data1$setName == "daccDietSelection_Rew_Delay"),]);
#   obs <- rbind(obs,fits$coefficients[2:3]);
#   ps <-rbind(ps,summary(fits)$coefficients[,4][2:3])
# }
# 
#obs <- NULL
#ps <- NULL
#obs2 <- NULL
#sigma <- NULL
#for (i in 1:125) {
#  obs2 <- rbind(obs2,runReg(data[which(data$neuron==i & data$setName == 'ofcSimulation_exp_desc'),]))
#}
# 
# q = runReg(data[ which(data$neuron=='1' 
#                        & data$setName == "ofcCuriosity_Rew_Info"), ])
