library(compiler)
library(iterators)
library(snow)
library(foreach)
library(doSNOW)
library(gbm)
source('functions.r')
library(dismo)
require(dismo)

nproc <- 3
cl.tmp = makeCluster(rep('localhost',nproc), type='SOCK')
registerDoSNOW(cl.tmp)
getDoParWorkers()

############################################################################
## Sensitivity analysis
############################################################################

#gravid lethal
sa_GL=read.csv("SArawresultsfinal_gravidlethal.csv", header=TRUE)

# convert breeding cycles to years
sa_GL$time           = (sa_GL$T_min-54)/6
# convert range to distance
sa_GL$D              = (sa_GL$d-1)/2
# initialize extinct
sa_GL$extinct=0

# sort out extincts before inoculation
for (i in 1:nrow(sa_GL)) {
  if( sa_GL[i,]$N_min == 0 & sa_GL[i,]$N_inoc > 0) { 
    sa_GL[i,]$extinct = 1
  } else if( sa_GL[i,]$N_inoc == 0){
    sa_GL[i,]$extinct = NA
  } else { 
    sa_GL[i,]$extinct = 0
  }
  print(i)
}

# subset extinct ones
sa_GL_ext = subset(sa_GL, extinct==1)

# subset by release strategies
sa_GL_BR = subset(sa_GL, strategy == "BD" & !is.na(extinct))
sa_GL_TCG = subset(sa_GL, strategy == "TC" & detection == "GD" & !is.na(extinct))
sa_GL_TCN = subset(sa_GL, strategy == "TC" & detection == "ND" & !is.na(extinct))

## number of simulations (nSims) and number of parameter samples (nSamples)
## note that this code assumes you are running a single simulation iteration per parameter sample (as recommended in the paper)
## hence nSims = nSamples
nSims <- nrow(sa_GL_ext)
nSamples <- nSims

########################################################################################
## Emulation
########################################################################################

## specify:
## the focus response variable - 'Extinct' (indicator or population persistence) or 'r' (population growth rate)
## tree.complexities to test with the BRT emulators (just testing no interactions [tc=1] and first-order interactions [tc=2] below)

tree.complexities <- 3:
  ## subsample sizes to test during the emulation step
  (subsamples <- seq(1000,nSims,by=1000))

#p extinction
fita_BR <- gbm.step(data=sa_GL_BR, gbm.x=4:10, gbm.y=17,family='bernoulli',tree.complexity=3,n.folds=5,tolerance.method='auto',max.trees=200000)
sensitivityplota = summary(fita_BR)
summary(fita_BR)

fita_TCG <- gbm.step(data=sa_GL_TCG, gbm.x=4:10, gbm.y=17,family='bernoulli',tree.complexity=3,n.folds=5,tolerance.method='auto',max.trees=200000)
sensitivityplota = summary(fita_TCG)
summary(fita_TCG)

fita_TCN <- gbm.step(data=sa_GL_TCN, gbm.x=4:10, gbm.y=17,family='bernoulli',tree.complexity=3,n.folds=5,tolerance.method='auto',max.trees=200000)
sensitivityplota = summary(fita_TCN)
summary(fita_TCN)



############################################################################
#Y_linked
sa_Y=read.csv("YeditBR-SAresultssep30.csv", header=TRUE)

# convert breeding cycles to years
sa_Y$time           = (sa_Y$T_min-54)/6
# convert range to distance
sa_Y$D              = (sa_Y$d-1)/2
# initialize extinct
sa_Y$extinct=0

# sort out extincts before inoculation
for (i in 1:nrow(sa_Y)) {
  if( sa_Y[i,]$N_min == 0 & sa_Y[i,]$N_inoc > 0) { 
    sa_Y[i,]$extinct = 1
  } else if( sa_Y[i,]$N_inoc == 0){
    sa_Y[i,]$extinct = NA
  } else { 
    sa_Y[i,]$extinct = 0
  }
  print(i)
}

# subset extinct ones
sa_Y_ext = subset(sa_Y, extinct==1)

# subset by release strategies
sa_Y_BR = subset(sa_Y, !is.na(extinct))

## number of simulations (nSims) and number of parameter samples (nSamples)
## note that this code assumes you are running a single simulation iteration per parameter sample (as recommended in the paper)
## hence nSims = nSamples
nSims <- nrow(sa_Y)
nSamples <- nSims

########################################################################################
## Emulation
########################################################################################

## specify:
## the focus response variable - 'Extinct' (indicator or population persistence) or 'r' (population growth rate)
## tree.complexities to test with the BRT emulators (just testing no interactions [tc=1] and first-order interactions [tc=2] below)

tree.complexities <- 3:
  ## subsample sizes to test during the emulation step
  (subsamples <- seq(1000,nSims,by=1000))

#p extinction
fit_Y_BR <- gbm.step(data=sa_Y_BR, gbm.x=4:10, gbm.y=17,family='bernoulli',tree.complexity=3,n.folds=5,tolerance.method='auto',max.trees=200000)
sensitivityplota = summary(fit_Y_BR)
summary(fit_Y_BR)



