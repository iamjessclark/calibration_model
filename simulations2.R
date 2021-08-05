require(truncnorm)
# simulating low prevalence G-score
modelDatalowgs <- function(nInd, prev, meanval, nFECrep=2){
  
  posts <- as.data.frame(cbind(c(kkGscore.model$mcmc[[1]][,"k"], kkGscore.model$mcmc[[2]][,"k"]),
                                      c(kkGscore.model$mcmc[[1]][,"intercept"], kkGscore.model$mcmc[[2]][,"intercept"]),
                                      c(kkGscore.model$mcmc[[1]][,"rtnb"], kkGscore.model$mcmc[[2]][,"rtnb"])))
  postsamples <- posts[sample.int(nrow(posts), 500, replace=FALSE, prob=NULL),]
  
  index <- runif(1,1,length(postsamples))
  k <- postsamples[index,1]
  inter <- postsamples[index,2]
  OvDisp <- postsamples[index,3]
  
  Status <- rbinom(nInd, size = 1, prob = prev)
  getoptim<-optimize(f=getgammalow,c(0, 1000000))
  #find the value of shape that is closest to satisfying a 95% quantile at 400
  #your params
  rate<-getoptim$minimum
  shape<-meanval*rate
  
  #Infectlvl <- rep(0, nInd)
  indexpos <- which(Status==1)
  #for(i in 1:length(indexpos)){
  #  Infectlvl[indexpos[i]] <- rgamma(1,shape = shape, rate = rate)
  #}
  
  Infectlvl <- rgamma(length(which(Status==1)),shape = shape, rate = rate)
  
  FECall <- matrix(rep(0,nInd*nFECrep),nrow=nFECrep)
  #lambdanums <- Infectlvl[which(Status==1)]
  #FECall[,which(Status==1)] <- sapply(lambdanums, 
                                      #function(x) {return(rnbinom(n=nFECrep,mu = x, size = OvDisp))})
  
  FECall[,which(Status==1)] <- sapply(Infectlvl, function(x) {return(rnbinom(n=nFECrep,mu = x, size = OvDisp))})
  
  # gscore
  indexneg <- which(Status==0)
  gscore <- rep(0,nInd)
  for(i in 1:length(indexneg)){
    gscore[indexneg[i]] <- round(rtruncnorm(n=1, mean=0, a=0, b=2, sd=sqrt(1/1.245443)),0)
  }
  
  for(i in 1:length(indexpos)){
    gscore[indexpos[i]] <- sapply((9 / (1 + exp(-k*(Infectlvl[i]-inter)))), 
                        function(x){round(rtruncnorm(n=1, mean=x, a=0, b=9, sd=sqrt(1/1.245443)),0)})
  }
  
  return(list(Status, FECall,Infectlvl, gscore))
}

# simulating low prevalence CCA
modelDatalowcca <- function(nInd, prev, meanval, nFECrep=2){
  
  posts <- as.data.frame(cbind(c(kkcca.model$mcmc[[1]][,"k"], kkcca.model$mcmc[[2]][,"k"]),
                               c(kkcca.model$mcmc[[1]][,"intercept"], kkcca.model$mcmc[[2]][,"intercept"]),
                               c(kkcca.model$mcmc[[1]][,"rtnb"], kkcca.model$mcmc[[2]][,"rtnb"])))
  postsamples <- posts[sample.int(nrow(posts), 500, replace=FALSE, prob=NULL),]
  
  index <- runif(1,1,length(postsamples))
  k <- postsamples[index,1]
  inter <- postsamples[index,2]
  OvDisp <- postsamples[index,3]
  
  Status <- rbinom(nInd, size = 1, prob = prev)
  getoptim<-optimize(f=getgammalow,c(0, 1000000))
  #find the value of shape that is closest to satisfying a 95% quantile at 400
  #your params
  rate<-getoptim$minimum
  shape<-meanval*rate
  
  #Infectlvl <- rep(0, nInd)
  indexpos <- which(Status==1)
  #for(i in 1:length(indexpos)){
  #  Infectlvl[indexpos[i]] <- rgamma(1,shape = shape, rate = rate)
  #}
  
  Infectlvl <- rgamma(length(which(Status==1)),shape = shape, rate = rate)
  
  FECall <- matrix(rep(0,nInd*nFECrep),nrow=nFECrep)
  #lambdanums <- Infectlvl[which(Status==1)]
  #FECall[,which(Status==1)] <- sapply(lambdanums, 
  #function(x) {return(rnbinom(n=nFECrep,mu = x, size = OvDisp))})
  
  FECall[,which(Status==1)] <- sapply(Infectlvl, function(x) {return(rnbinom(n=nFECrep,mu = x, size = OvDisp))})
  
  #cca
  CCA <- rep(0,nInd)
  indexneg <- which(Status==0)
  for(i in 1:length(indexneg)){
    #CCAall[i] <- round(rtruncnorm(n=1, mean=0, a=0, b=1, sd=sqrt(1/3.093451)),0)
    CCA[indexneg[i]] <- round(rtruncnorm(n=1, mean=0, a=0, b=1, sd=sqrt(1/3.47182)),0)
  }
  
  for(i in 1:length(indexpos)){
    CCA[indexpos[i]] <- sapply((4 / (1 + exp(-k*(Infectlvl[i]-inter)))), 
                               function(x){round(rtruncnorm(n=1, mean=x, a=0, b=4, sd=sqrt(1/3.47182)),0)})
  }
  
  return(list(Status, FECall,Infectlvl, CCA))
}

# running sim repeats

runRepeatsgs <- function(index){
  return(apply(replicate(n=nRep,modelDatalowgs(nInd, prev=prev, meanval=meanval, nFECrep = nFECrep)),2,as.list))
}

runRepeatscca <- function(index){
  return(apply(replicate(n=nRep,modelDatalowcca(nInd, prev=prev, meanval=meanval, nFECrep = nFECrep)),2,as.list))
}

# optim function
getgammalow<-function(rate){
  #this function takes a shape value and your given mean value and works out the 95% quantile. It then returns the absolute difference between that and 400.
  
  shape<-meanval*rate
  
  
  return(abs((qgamma(p=c(0.99),shape=shape,rate=rate)-16.7)))
}


# extracting prevalences from lists 

extractprevs <- function(input, prevalences, dim, threshold, nInd, nRuns){
  prevs <- list()
  
  for(i in 1:length(input)){
    prevs[[i]] <- list()
  }
  
  for(i in 1:length(prevalences)){
    for(j in 1:nRuns){
    prevs[[i]][[j]] <- length(which(input[[i]][[j]][[dim]]>=threshold))/nInd
    }
  }
  prevs <- as.vector(do.call(rbind,(flatten((prevs)))))
  return(prevs)
}


extractkkprevs <- function(input, prevalences, dim, threshold, nInd, nRuns){
  prevs <- list()
  
  for(i in 1:length(input)){
    prevs[[i]] <- list()
  }
  
  for(i in 1:length(prevalences)){
    for(j in 1:nRuns){
      prevs[[i]][[j]] <- length(which(colSums(as.data.frame(input[[i]][[j]][[2]]))/2>0))/nInd
    }
  }
  
  prevs <- as.vector(do.call(rbind,(flatten((prevs)))))
  return(prevs)
}

# extracting kkcounts and scores

getVals <- function(input){
  
  kk <- list()
  cca <- list()
  index <- sample(1:100, 50, replace=F)
  for(i in 1:50){
      kk[[i]] <- colSums(as.data.frame(input[[index[i]]][[index[i]]][[2]]))/2
      cca[[i]] <- unlist(input[[index[i]]][[index[i]]][[4]])
    }
 
 kk <- t(do.call(rbind, kk) )
 cca <- t(do.call(rbind, cca))
 output <- list(kk, cca)
 return(output)
}

