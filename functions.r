##############################################################################################
## 02. Functions.R
## This script specifies the following functions:
## (1) SA.func - runs a global sensitivity analysis for population model
## (2) emulation.func - emulates the sensitivity analysis output with boosted regression trees
##############################################################################################



##########################################################################################
## Sensitivity Analysis Function
## This function requires the package 'lhs', and specification of the following arguments:
## SAvars - vector of parameter names to be included
## nSims - total number of simulations
## nSamples - number of parameter samples to draw
## type - type of sampling to implement, currently must be 'random' or 'latin'
## NOTE: the number of model iterations per sample (nIter) is calculated as nSims/nSample
##########################################################################################

SA.func <- function(SAvars, nSims, nSamples, type, ...) {
  
  dir.create('sensitivity_analysis', showWarnings = F)
  
  ## load  packages
  if(type=='latin') library(lhs)
  
  ## calculate number of iterations per sample
  nIter <- nSims/nSamples
  
  ####################################################
  ## Generate parameter samples
  ####################################################
  
  ## first set up template with default values
  samples <- expand.grid(parDefault.list)
  samples <- samples[rep(seq_len(nrow(samples)), each=nSamples),]
  
  nVars <- length(SAvars)
  
  ## generate uniform samples between 0 and 1 for required parameters
  if (type=='random') {
    raw.samples <- matrix(NA, nrow=nSamples, ncol=nVars)
    for(i in 1:nVars) {raw.samples[,i] <- runif(nSamples,0,1)}
  } else if (type=='latin') {
    raw.samples <- randomLHS(n=nSamples, k=nVars)
  }
  
  ## transform using required ranges
  for(i in 1:length(SAvars)) {
    if (SAvars[i]=='nStart') {
      temp.nStart.vec <- qunif(raw.samples[,i], min=parRange.list[[SAvars[i]]][1], max=parRange.list[[SAvars[i]]][2])
      samples[,SAvars[i]] <- (temp.nStart.vec%/%2)*2 + ifelse(temp.nStart.vec%%2<1,0,2)
    } else if (SAvars[i]%in%c('ageMaturity','P')) {
      samples[,SAvars[i]] <- round(qunif(raw.samples[,i], min=parRange.list[[SAvars[i]]][1]-0.5, max=parRange.list[[SAvars[i]]][2])+0.5)
    } else {
      samples[,SAvars[i]] <- qunif(raw.samples[,i], min=parRange.list[[SAvars[i]]][1], max=parRange.list[[SAvars[i]]][2])
    }
  }
  
  samples$nSims <- nSims
  samples$nSamples <- nSamples
  samples$nIter <- nIter
  
  ## run the population projection
  result <- foreach(row=1:nrow(samples), .combine=rbind, .packages=c('MASS'), .export=c('popProj.func.comp')) %dopar% {
    popProj.func.comp(nStart=samples$nStart[row],ageMaturity=samples$ageMaturity[row],sr=samples$sr[row],m=samples$m[row],    
                      m.Imult=samples$m.Imult[row],s0.mult=samples$s0.mult[row],s1plus=samples$s1plus[row],
                      s0.Imult=samples$s0.Imult[row],s1plus.Imult=samples$s1plus.Imult[row],pBreed=samples$pBreed[row],
                      pBreed.Imult=samples$pBreed.Imult[row],pMatTrans=samples$pMatTrans[row],pIoutside=samples$pIoutside[row],
                      beta=samples$beta[row],pRecover=samples$pRecover[row],pResistant=samples$pResistant[row],pLoseResistance=samples$pLoseResistance[row],
                      aPred=samples$aPred[row],hPred=samples$hPred[row],P=samples$P[row],nIter=samples$nIter[row]) 
  }
  
  ## save the results
  save.nm <- paste('sensitivity_analysis/nSims=',nSims,'_nSamples=',nSamples,'_nIter=',nIter,'_type=',type,'.results',sep='')
  assign(save.nm, result)
  save(list=save.nm, file=save.nm)
  
  ## return the results
  return(result)
}

################################################################################################################################
## Emulation Function
## This function emulates the sensitivity analysis output with boosted regression trees with different interaction depths
## It takes the following arguments:
## data - the sensitivity analysis output to use (produced by the function SA.func above)
## SAvars - vector of parameter names to be included
## resp - the focal response variable
## subsample - vector of subsamples (i.e., number of data rows) for which emulation will be performed
## tree.complexities - vector of tree complexities (interaction depths) to test
################################################################################################################################

emulation.func <- function(data, SAvars, resp, subsample, tree.complexities, ...) {
  
  require(dismo)
  
  ## clean up and create folders
  gc()
  dir.create('emulation',showWarnings=F)
  dir.create('emulation/brts',showWarnings=F)
  dir.create('emulation/results',showWarnings=F)
  
  ## subset required simulation results
  dataset <- data.frame(data[1:subsample,])
  # dataset$r[dataset$r%in%c('-Inf','Inf')] <- NA
  
  ## statistical distribution for fitting BRTs
  
  brt.dist <- ifelse(resp=='Extinct','bernoulli','gaussian')
  ##brt.dist <- ifelse(resp=='cycle','poisson')
  
  ## fit BRT emulators of different tree complexities for given response variable
  for (i in 1:length(tree.complexities)) {
    
    tc <- tree.complexities[i]
    
    ## BRT 
    brt.fit <- NULL
    if((resp=='Extinct' & length(unique(dataset$Extinct))>1) | resp=='r') {
    ##if((resp=='cycle' & length(unique(dataset$cycle))>1) | resp=='r') {
      x.col <- which(names(dataset)%in%SAvars)
      y.col <- which(names(dataset)==resp)
      brt.dataset <- dataset[!is.na(dataset[,y.col]),]
      brt.fit <- try(gbm.step(data=brt.dataset, gbm.x=x.col, gbm.y=y.col,family=brt.dist,tree.complexity=i,n.folds=5,tolerance.method='auto',max.trees=200000))
      ## try again with a decreased learning rate if necessary
      indic <- 1
      while((indic <= 49) & (inherits(brt.fit,'try-error')|is.null(brt.fit))) {
        lr <- 0.01*0.9
        brt.fit <- try(gbm.step(data=brt.dataset, gbm.x=x.col, gbm.y=y.col,family=brt.dist,tree.complexity=i,n.folds=5,tolerance.method='auto',tolerance=tol,learning.rate=lr,step.size=ss,max.trees=200000))
        indic <- indic + 1
      } 
      if (!inherits(brt.fit,'try-error') & !is.null(brt.fit)) {
        ## save brt file
        (brt_save.nm <- paste('subsample.',subsample,'.brt',tc,'.',resp,sep=''))
        assign(brt_save.nm,brt.fit)
        save(list=brt_save.nm,file=paste('emulation/brts/',brt_save.nm,sep=''))
        rm(list=c(brt_save.nm,'brt.fit'))
      }
      gc()
      
    }
  }
  return('Done')
}

###################################################################################################################################
## Emulation Summary Function
## This function summarises the emulation results produced by 'emulation.func' above
## It returns a list of 2 data frames that store:
## the cross-validation deviance for each emulation
## the stability of relative influence metrics results as the subsample size is increased
###################################################################################################################################

emulation.summary.func <- function(resp=resp, subsamples=subsamples, tree.complexities=tree.complexities) {
  
  require(MDM)

  cvDev.df <- betaDiv.df <- data.frame(matrix(NA, nrow=length(subsamples), ncol=length(tree.complexities)+1))
  names(cvDev.df) <- names(betaDiv.df) <- c('subsample',paste('tc',tree.complexities,sep='.'))
  cvDev.df$subsample <- betaDiv.df$subsample <- subsamples
  ri.df <- NULL
  
  for (i in 1:length(tree.complexities)) {
    
    tc <- tree.complexities[i]
    
    file.nms <- paste('subsample.',subsamples,'.brt',tc,'.',resp,sep='')
    ri.temp <- data.frame(Parameter=sort(SAvars))
    
    for (j in 1:length(file.nms)) {
      
      load(paste('emulation/brts/',file.nms[j],sep=''))
      
      ## Relative influence part
      brt <- eval(as.name(file.nms[j]))
      ri <- brt$contributions
      ri <- ri[order(ri$var),]
      ri.temp <- cbind(ri.temp,ri[,2,drop=F])
      
      ## cv deviance part
      dev <- brt$cv.statistics$deviance.mean
      cvDev.df[j,paste('tc',tc,sep='.')] <- dev
      
      ## Clean up
      rm(list=c(file.nms[j],'brt'))
      print(j)
    }
    
    ri.temp$tc <- tc
    rownames(ri.temp) <- 1:nrow(ri.temp)
    ri.df <- rbind(ri.df,ri.temp)
  }
  
  ## Calculate beta diversities between pairs
  for (i in 1:length(tree.complexities)) {
    tc <- tree.complexities[i]
    df <- ri.df[ri.df$tc==tc,]
    betaDiv.dummy <- data.frame(subsample=subsamples,resp=NA)
    for (j in 2:length(subsamples)) {
      beta.div <- as.numeric(ed(t(as.matrix(df[,j:(j+1)])),q=1,retq=T)['beta'])
      betaDiv.dummy[j,2] <- beta.div    
    }
    betaDiv.df[,paste('tc',tc,sep='.')] <- betaDiv.dummy[,2]
  }

  return(list(cvDev = cvDev.df, betaDiv=betaDiv.df))
}

###################################################################################################################################
## Emulation Plotting Function
## This function creates a pdf plot of the results produced the 'emulation summary'
## the plot.name argument allows you to name the pdf
###################################################################################################################################

emulationPlot.func <- function(data, plot.name) {
  
  cvDev <- data[['cvDev']]
  betaDiv <- data[['betaDiv']]
  
  pdf(paste(plot.name,'pdf',sep='.'),width=6,height=8)
  par(mfrow=c(2,1),mar=c(4,5,2,2),oma=c(1,1,1,1),mgp=c(2.4,0.6,0))
  ylimits <- range(cvDev[,2:ncol(cvDev)],na.rm=T)
  plot(cvDev[,1:2],type='o',ylab='Cross-validation deviance',xlab='',pch=21,bg=1,lty=2,ylim=ylimits)
  for (i in 3:ncol(cvDev)) lines(cvDev[,c(1,i)],col=i-1,type='o',pch=21,bg=i-1,lty=2)
  legend('topright',gsub('tc.','',names(cvDev)[-1]),col=1:(ncol(cvDev)-1),lty=1,pch=21,pt.bg=1:(ncol(cvDev)-1), title='Tree Complexity')
  ylimits <- range(betaDiv[,2:ncol(betaDiv)],na.rm=T)
  plot(betaDiv[,1:2],type='o',ylab='Beta-diversity of\nrelative influence metrics',xlab='Subsample Size',pch=21,bg=1,lty=2,ylim=ylimits)
  for (i in 3:ncol(betaDiv)) lines(betaDiv[,c(1,i)],col=i-1,type='o',pch=21,bg=i-1,lty=2)
  legend('topright',gsub('tc.','',names(betaDiv)[-1]),col=1:(ncol(betaDiv)-1),lty=1,pch=21,pt.bg=1:(ncol(betaDiv)-1), title='Tree Complexity')
  dev.off()
  
}

