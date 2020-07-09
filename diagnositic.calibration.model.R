# Calibrating Schistosomiasis Diagnostics
# Jessica Clark
# Using the bayesian framework developed in the HMM to look at each time point independently
# Status will not be tracked across time points but will be estimated independently

# this is for just the KK data
CalModKK <- function (N, Ti, R, KK) {
  ## Set seed ##
  .RNG.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  
  ## Initialize Status ##
  Status <- matrix(1,nrow=210, ncol=4)
  sigma <-  0.001
  
  m <- "model {
    
    # Prior prevalence  #

    prev ~ dbeta(1,1)
    
    # Prior KKs #
    
    rtnb ~ dgamma(286.3235, 223.4650)
    sigma ~ dgamma(0.001, 0.001)
    
    tkksh ~ dgamma(83.88592,125.70331)
    tkkrt1 ~ dbeta(47.13542,633.08366)
    tkkrt <- tkkrt1/(1-tkkrt1)
    
    # model 
    for(n in 1:N){	# run through pop
  
      for (t in 1:Ti){ # run through time
        Status[n,t] ~ dbern(prev)
        
        tKK[n,1,t] <- 0
        tKK[n,2,t] ~ dgamma(tkksh, tkkrt)
        
        lambda[n,t] <- tKK[n,Status[n,t]+1,t] 
          
          for(r in 1:R){ # run through repeated measures to set the baseline KK over the repeated measures, this is the LLH
             KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb)
          } # end or r loop
        } # end time
      } # end n loop
    
    #inits# .RNG.seed, .RNG.name, Status, sigma
    #data# N, Ti, R, KK
    #monitor# prev, rtnb, Status
}"
  
  
  # Run model #
  Results <- run.jags(m, burnin=10000, sample=20000, n.chains=2, jags.refresh = 1, method = 'parallel',
                      plots = F, silent.jags = F)
  return(Results)
}

#### CCA and KK ####

CalModKKCCA <- function (N, Ti, R, KK, CCA) {
  ## Set seed ##
  ##this is so that they start on different values ##
  .RNG.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  
  ## Initialize Status ##
  Status <- matrix(1,nrow=210, ncol=4)
  sigma <-  0.001
  prob <- array(rep(CCA,2),dim=c(N,Ti,2))
  
  m <- "model {
  # Prior prevalence / clearance / reinfection #
  prev ~ dbeta(1,1)

  # Prior KKs 
  rtnb ~ dgamma(286.3235, 223.4650)
      
  tkksh ~ dgamma(83.88592,125.70331)
  tkkrt1 ~ dbeta(47.13542,633.08366)
  tkkrt <- tkkrt1/(1-tkkrt1)

  k ~ dgamma(0.001, 0.001)
  intercept ~ dgamma(0.001, 0.001)
  
  # MODEL 
  for(n in 1:N){	# run through pop
    
    for (t in 1:Ti){ # run through timepoints
      
      Status[n,t] ~ dbern(prev)
      
      lambda[n,t] <- tKK[n,Status[n,t]+1,t] # this needs to go here because each time the r loop reiterates it will replace this value
        
        for( r in 1:R){  # run through repeat measures
          KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb) # generating the data with noise and then sampling from the dataset with the gamma sh/rt parameters?
      } # end or r loop

    
      tKK[n,1,t] <- 0                                                                                                                                                                                                                                                                                               
      tKK[n,2,t] ~ dgamma(tkksh, tkkrt)
      
      prob[n,t,1] ~ dnorm(0,3.093451)
      prob[n,t,2] ~ dnorm(4 / (1 + exp(-k*(tKK[n,2,t]-intercept))),3.093451)
       
      CCA[n,t] ~ dround(prob[n,t,Status[n,t]+1],0)
      
      } # end timestep 1 t loop
    } #N
  
  #inits# .RNG.seed, .RNG.name, Status, sigma, prob   
  #data# N, Ti, R, KK, CCA
  #monitor#  rtnb, prev, k, intercept, Status
}"
  
  # Run model #
  Results <- run.jags(m, burnin=10000, sample=20000, n.chains=2, jags.refresh = 1, method = 'parallel',
                      plots = F, silent.jags = F)
  return(Results)
}

#### KK and G-Score model ####

CalModKKGScore <- function (N, Ti, R, KK, CCA10) {
  ## Set seed ##
  ##this is so that they start on different values ##
  .RNG.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  
  ## Initialize Status ##
  Status <- matrix(1,nrow=210, ncol=4)
  prob <- array(rep(CCA10,2),dim=c(N,Ti,2))
  sigma<- 0.001
  
  m <- "model {
  
  # MODEL 
   # Prior prevalence #
    prev ~ dbeta(1,1)
   
    # Prior KKs 
    rtnb ~ dgamma(286.3235, 223.4650)

    tkksh ~ dgamma(0.001, 0.001)
    tkkrt1 ~ dbeta(1,1)
    tkkrt <- tkkrt1/(1-tkkrt1)
    
    k ~ dgamma(0.001, 0.001)
    intercept ~ dgamma(0.001, 0.001)
    
    for(n in 1:N){	# run through pop
      
      for (t in 1:Ti){ # run through timepoints
      
            Status[n,t] ~ dbern(prev)
            
            lambda[n,t] <- tKK[n,Status[n,t]+1,t] 
            
            for( r in 1:R){  # run through repeat measures
              KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb) 
              } # end or r loop
      
            tKK[n,1,t] <- 0
            tKK[n,2,t] ~ dgamma(tkksh, tkkrt)
  
            prob[n,t,1] ~ dnorm(0,1.093606) 
            prob[n,t,2] ~ dnorm(9 / (1 + exp(-k*(tKK[n,2,t]-intercept))),1.093606)
            
            CCA10[n,t] ~ dround(prob[n,t,Status[n,t]+1],0)
            
      } # t loop
    } # N loop   
    
 
  #inits# .RNG.seed, .RNG.name, Status, prob,  sigma
  #data# N, Ti, R, KK, CCA10
  #monitor#  prev, rtnb, k, intercept, Status
  
}"
  
  # Run model #
  Results <- run.jags(m, burnin=10000, sample=20000, n.chains=2, jags.refresh = 1, method = 'parallel',
                      plots = F, silent.jags = F)
  return(Results)
}



#### Data manipulation functions ####


#### Data loading ####

#### Kato Katz Data ####

loadKKdata <- function(nameFile){
  dt <- read.csv(nameFile) # read data file
  
  # name different weeks
  dt$dateN <- NA
  dt$dateN[which(dt$date=="25/09/2017" | dt$date=="26/09/2017" | dt$date=="27/09/2017"
                 | dt$date=="28/09/2017" | dt$date=="29/09/2017" | dt$date=="02/10/2017")] <- "Pre-T"
  
  dt$dateN[which(dt$date=="23/10/2017" | dt$date=="24/10/2017" | dt$date=="25/10/2017"
                 | dt$date=="26/10/2017" | dt$date=="27/10/2017")] <- "3 weeks"
  
  dt$dateN[which(dt$date=="04/12/2017" | dt$date=="05/12/2017")] <- "9 weeks"
  
  dt$dateN[which(dt$date=="01/03/2018" | dt$date=="05/03/2018" | dt$date=="06/03/2018"
                 | dt$date=="07/03/2018" | dt$date=="08/03/2018" | dt$date=="09/03/2018")] <- "6 months"
  
  dates <- c(na.omit(unique(dt$dateN)))
  # set as array dimension (number children, number timesteps, number max of repeats)
  children <- unique(dt$child_id)
  KK <- array(NA,dim = c(length(children),4,6))
  for (i in 1:length(children)){
    KK[i,,] <- getKKChild(children[i],dt,dates)
  }
  return(KK)
}

getKKChild <- function(ID,dt,dates){
  dts <- subset(dt,child_id==ID,select = c(child_id,Sm_A,Sm_B,dateN))
  KKChild <- matrix(NA,nrow=length(dates),ncol=6)
  for (i in 1:length(dates)){
    KKChild[i,] <- getKKChildWeek(dates[i],dts)
  }
  return(KKChild)
}

getKKChildWeek <- function(weekN,dts){
  repeatsKK <- c(dts$Sm_A[which(dts$dateN==weekN)],dts$Sm_B[which(dts$dateN==weekN)])
  length(repeatsKK)<-6 #set max number of repeats!
  return(repeatsKK)
}

getKKChildIDs <- function(nameFile){
  dt <- read.csv(nameFile)
  return(sort(unique(dt$child_id)))
}

#### Normal CCA Data ####

loadCCAdata <- function(nameFile){
  dt <- read.csv(nameFile, stringsAsFactors = F)
  dt$Poppy[which(dt$Poppy == "I*" | dt$Poppy=="-")]<-NA
  
  # restructure for categorical
  dt$Poppy[which(dt$Poppy=="T")] <- 1
  dt$Poppy[which(dt$Poppy==3)]<-4
  dt$Poppy[which(dt$Poppy==2)]<-3
  dt$Poppy[which(dt$Poppy==1)]<-2
  dt$Poppy[which(dt$Poppy==0.5)]<-1
  
  # name different weeks
  dt$dateN <- NA
  dt$dateN[which(dt$date_on_tube=="25/09/2017" | dt$date_on_tube=="26/09/2017" | dt$date_on_tube=="27/09/2017"
                 | dt$date_on_tube=="28/09/2017" | dt$date_on_tube=="29/09/2017" | dt$date_on_tube=="02/10/2017")] <- "Pre-T"
  
  dt$dateN[which(dt$date_on_tube=="23/10/2017" | dt$date_on_tube=="24/10/2017" | dt$date_on_tube=="25/10/2017"
                 | dt$date_on_tube=="26/10/2017" | dt$date_on_tube=="27/10/2017")] <- "3 weeks"
  
  dt$dateN[which(dt$date_on_tube=="04/12/2017" | dt$date_on_tube=="05/12/2017")] <- "9 weeks"
  
  dt$dateN[which(dt$date_on_tube=="01/03/2018" | dt$date_on_tube=="05/03/2018" | dt$date_on_tube=="06/03/2018"
                 | dt$date_on_tube=="07/03/2018" | dt$date_on_tube=="08/03/2018" | dt$date_on_tube=="09/03/2018")] <- "6 months"
  
  dates <- c(na.omit(unique(dt$dateN)))
  # set as array dimension (number children, number timesteps)
  children <- unique(dt$CID)
  CCA <- array(NA,dim = c(length(children),4))
  for (i in 1:length(children)){
    CCA[i,] <- getCCAChild(children[i],dt,dates)
  }
  return(CCA)
} # end of  function 

# function for getting CCA child
getCCAChild <- function(ID,dt,dates){
  dts <- subset(dt,CID==ID,select = c(CID,Poppy,dateN))
  CCAChild <- rep(NA,length(dates))
  for (i in 1:length(dates)){
    if (!is.null(which(dts$dateN==dates[i]))){
      CCAChild[i] <- ifelse(is.null(dts$Poppy[which(dts$dateN==dates[i])]),NA,as.numeric(dts$Poppy[which(dts$dateN==dates[i])]))
    }
  }
  return(CCAChild)
}

# function for getting CCA child ID
getCCAChildIDs <- function(nameFile){
  dt <- read.csv(nameFile)
  return(sort(unique(dt$CID)))
}

#### GScore data ####

loadgscoredata <- function(nameFile){
  #restructure so that it is the right scoring so that a
  # score of 1 which is negative, is given a score of 0.
  dt <- read.csv(nameFile, stringsAsFactors = F)
  
  dt$gscore <- dt$gscore-1
  
  # name different weeks
  dt$dateN <- NA
  dt$dateN[which(dt$date_on_tube=="25/09/2017" | dt$date_on_tube=="26/09/2017" | dt$date_on_tube=="27/09/2017"
                 | dt$date_on_tube=="28/09/2017" | dt$date_on_tube=="29/09/2017" | dt$date_on_tube=="02/10/2017")] <- "Pre-T"
  
  dt$dateN[which(dt$date_on_tube=="23/10/2017" | dt$date_on_tube=="24/10/2017" | dt$date_on_tube=="25/10/2017"
                 | dt$date_on_tube=="26/10/2017" | dt$date_on_tube=="27/10/2017")] <- "3 weeks"
  
  dt$dateN[which(dt$date_on_tube=="04/12/2017" | dt$date_on_tube=="05/12/2017")] <- "9 weeks"
  
  dt$dateN[which(dt$date_on_tube=="01/03/2018" | dt$date_on_tube=="05/03/2018" | dt$date_on_tube=="06/03/2018"
                 | dt$date_on_tube=="07/03/2018" | dt$date_on_tube=="08/03/2018" | dt$date_on_tube=="09/03/2018")] <- "6 months"
  
  dates <- c(na.omit(unique(dt$dateN)))
  # set as array dimension (number children, number timesteps)
  children <- unique(dt$CID)
  CCAgscore <- array(NA,dim = c(length(children),4))
  for (i in 1:length(children)){
    CCAgscore[i,] <- getCCAChildgscore(children[i],dt,dates)
  }
  return(CCAgscore)
} # end of  function 

# function for getting CCA child
getCCAChildgscore <- function(ID,dt,dates){
  dts <- subset(dt,CID==ID,select = c(CID,gscore,dateN))
  CCAChildgscore <- rep(NA,length(dates))
  for (i in 1:length(dates)){
    if (!is.null(which(dts$dateN==dates[i]))){
      CCAChildgscore[i] <- ifelse(is.null(dts$gscore[which(dts$dateN==dates[i])]),NA,as.numeric(dts$gscore[which(dts$dateN==dates[i])]))
    }
  }
  return(CCAChildgscore)
}

# function for getting CCA child ID
getCCAChildIDsgscore <- function(nameFile){
  dt <- read.csv(nameFile)
  return(sort(unique(dt$CID)))
}

#### sort model output into time steps ####

time.steps <- function(model.output){
  require(reshape2)
  t1 <- as.data.frame(model.output[,1:210])
  t2 <- as.data.frame(model.output[,211:420])
  t3 <- as.data.frame(model.output[,421:630])
  t4 <- as.data.frame(model.output[,631:840])
  
  t1[nrow(t1)+1,] <- colSums(t1)
  t1[nrow(t1),] <- t1[nrow(t1),]/20000
  
  t1.means <- as.data.frame(t1[20001,1:210])
  
  t1.means$time <-as.factor("Baseline")
  
  T1 <- melt(t1.means)
  
  t2[nrow(t2)+1,] <- colSums(t2)
  t2[nrow(t2),] <- t2[nrow(t2),]/20000
  
  t2.means <- as.data.frame(t2[20001,1:210])
  
  t2.means$time <-as.factor("ThreeWeeks")
  
  T2 <- melt(t2.means)
  
  t3[nrow(t3)+1,] <- colSums(t3)
  t3[nrow(t3),] <- t3[nrow(t3),]/20000
  
  t3.means <- as.data.frame(t3[20001,1:210])
  
  t3.means$time <- as.factor("NineWeeks")
  
  T3 <- melt(t3.means)
  
  t4[nrow(t4)+1,] <- colSums(t4)
  t4[nrow(t4),] <- t4[nrow(t4),]/20000
  
  t4.means <- as.data.frame(t4[20001,1:210])
  
  t4.means$time <- as.factor("SixMonths")
  
  T4 <- melt(t4.means)
  
  props <- rbind(T1, T2, T3, T4)
  
  props$time <- as.factor(props$time)
  props$variable <- as.factor(props$variable)
  
  return(props)
}

# for getting the mean kk estimates out of the model output

time.pointmeans <- function(list.name){
  for(i in 1:length(list.name)){
    colnames(list.name[[i]]) <- CID
  }
  
  
  for(i in 1:length(list.name)){
    list.name[[i]] <- list.name[[i]] %>% pivot_longer(everything(),
                                                      names_to = "CID", values_to = "count")
  }
  
  counts <- cbind(list.name[[1]], list.name[[2]][,2], list.name[[3]][,2], list.name[[4]][,2], list.name[[5]][,2], list.name[[6]][,2]) 
  colnames(counts) <- c("child", "count1", "count2", "count3", "count4", "count5", "count6")
  counts$mean <- rowMeans(counts[,2:ncol(counts)])
  
  return(counts)
  
}

#### log curve function ####

logcurve <- function(x,k,intercept){
  y <- 1 / (1 + exp(-k*(x-intercept)))
  return(y)
}

#### add legend function ####

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

#### roc curve function ####

diag_spec <- function(status.output, mergekkcca, data.condition.column ){
  
  t1 <- as.data.frame(status.output[,1:210])
  t2 <- as.data.frame(status.output[,211:420])
  t3 <- as.data.frame(status.output[,421:630])
  t4 <- as.data.frame(status.output[,631:840])
  
  colnames(t1) <- CID
  colnames(t2) <- CID
  colnames(t3) <- CID
  colnames(t4) <- CID
  
  t1status.sample <- t1[sample.int(nrow(t1), 500, replace=FALSE, prob=NULL),]
  t1status.sample <- t(t1status.sample)
  t1status.sample <- as.data.frame(cbind(rownames(t1status.sample), data.frame(t1status.sample, row.names=NULL)))
  t1status.sample <- t1status.sample %>% rename(CID=`rownames(t1status.sample)`)
  
  # join to the kkpos
  merge <- mergekkcca %>% select(dateN, CID, data.condition.column) %>%
    filter(dateN=="Pre-T")%>%
    right_join(t1status.sample, by="CID")
  
  merge <- merge[-which(is.na(merge$dateN)==T),]
  merge <- merge[!duplicated(merge[,"CID"]),]  
  
  # get the sensitivity/ specificity 
  for(i in 1:nrow(merge)){
    for(c in 4:ncol(merge)){
      if(merge[i,3]==1 & merge[i,c]==1){
        merge[i,c] <- "TP"
      } else if(merge[i,3]==1 & merge[i,c]==0){
        merge[i,c] <- "FP"
      } else if(merge[i,3]==0 & merge[i,c]==0){
        merge[i,c] <- "TN" 
      } else {
        merge[i,c] <- "FN" 
      }
    }
  }
  merge <- merge %>% mutate_if(is.character, as.factor)
  speclist <- list()
  
  # calculate how often each one occurs and get proportions 
  for(i in 4:ncol(merge)){
    speclist[[i-3]] <-  merge %>% group_by(merge[,i]) %>% tally() %>%
      mutate(proportion = n/sum(n))
  }
  
  test.specs<- list()
  for(i in 1:length(speclist)){
    test.specs[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
    colnames(test.specs[[i]]) <- c("TPR", "FPR")
    test.specs[[i]][1,1] <- speclist[[i]][4,2]/(speclist[[i]][4,2]+speclist[[i]][1,2])
    test.specs[[i]][1,2] <- speclist[[i]][2,2]/(speclist[[i]][3,2]+speclist[[i]][2,2])
  }
  
  roct1 <- rbindlist(test.specs)
  roct1 <- roct1 %>% mutate(time = "pre-T")%>%mutate_if(is.character, as.factor)
  
  # 2nd time point 
  
  t2status.sample <- t2[sample.int(nrow(t2), 500, replace=FALSE, prob=NULL),]
  t2status.sample <- t(t2status.sample)
  t2status.sample <- as.data.frame(cbind(rownames(t2status.sample), data.frame(t2status.sample, row.names=NULL)))
  t2status.sample <- t2status.sample %>% rename(CID=`rownames(t2status.sample)`)
  
  
  merge <- mergekkcca %>% select(dateN, CID, data.condition.column) %>%
    filter(dateN=="3 weeks")%>%
    right_join(t2status.sample, by="CID")
  
  merge <- merge[-which(is.na(merge$dateN)==T),]
  merge <- merge[!duplicated(merge[,"CID"]),]  
  
  # get the sensitivity/ specificity 
  for(i in 1:nrow(merge)){
    for(c in 4:ncol(merge)){
      if(merge[i,3]==1 & merge[i,c]==1){
        merge[i,c] <- "TP"
      } else if(merge[i,3]==1 & merge[i,c]==0){
        merge[i,c] <- "FP"
      } else if(merge[i,3]==0 & merge[i,c]==0){
        merge[i,c] <- "TN" 
      } else {
        merge[i,c] <- "FN" 
      }
    }
  }
  merge <- merge %>% mutate_if(is.character, as.factor)
  speclist <- list()
  
  # calculate how often each one occurs and get proportions 
  for(i in 4:ncol(merge)){
    speclist[[i-3]] <-  merge %>% group_by(merge[,i]) %>% tally() %>%
      mutate(proportion = n/sum(n))
  }
  
  test.specs<- list()
  for(i in 1:length(speclist)){
    test.specs[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
    colnames(test.specs[[i]]) <- c("TPR", "FPR")
    test.specs[[i]][1,1] <- speclist[[i]][4,2]/(speclist[[i]][4,2]+speclist[[i]][1,2])
    test.specs[[i]][1,2] <- speclist[[i]][2,2]/(speclist[[i]][3,2]+speclist[[i]][2,2])
  }
  
  roct2 <- rbindlist(test.specs)
  roct2 <- roct2 %>% mutate(time = "3 weeks")%>%mutate_if(is.character, as.factor)
  
  # 3rd time point 
  
  t3status.sample <- t3[sample.int(nrow(t3), 500, replace=FALSE, prob=NULL),]
  t3status.sample <- t(t3status.sample)
  t3status.sample <- as.data.frame(cbind(rownames(t3status.sample), data.frame(t3status.sample, row.names=NULL)))
  t3status.sample <- t3status.sample %>% rename(CID=`rownames(t3status.sample)`)
  
  
  merge <- mergekkcca %>% select(dateN, CID, data.condition.column) %>%
    filter(dateN=="9 weeks")%>%
    right_join(t3status.sample, by="CID")
  
  merge <- merge[-which(is.na(merge$dateN)==T),]
  merge <- merge[!duplicated(merge[,"CID"]),]  
  
  # get the sensitivity/ specificity 
  for(i in 1:nrow(merge)){
    for(c in 4:ncol(merge)){
      if(merge[i,3]==1 & merge[i,c]==1){
        merge[i,c] <- "TP"
      } else if(merge[i,3]==1 & merge[i,c]==0){
        merge[i,c] <- "FP"
      } else if(merge[i,3]==0 & merge[i,c]==0){
        merge[i,c] <- "TN" 
      } else {
        merge[i,c] <- "FN" 
      }
    }
  }
  merge <- merge %>% mutate_if(is.character, as.factor)
  speclist <- list()
  
  # calculate how often each one occurs and get proportions 
  for(i in 4:ncol(merge)){
    speclist[[i-3]] <-  merge %>% group_by(merge[,i]) %>% tally() %>%
      mutate(proportion = n/sum(n))
  }
  
  test.specs<- list()
  for(i in 1:length(speclist)){
    test.specs[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
    colnames(test.specs[[i]]) <- c("TPR", "FPR")
    test.specs[[i]][1,1] <- speclist[[i]][4,2]/(speclist[[i]][4,2]+speclist[[i]][1,2])
    test.specs[[i]][1,2] <- speclist[[i]][2,2]/(speclist[[i]][3,2]+speclist[[i]][2,2])
  }
  
  roct3 <- rbindlist(test.specs)
  roct3 <- roct3 %>% mutate(time = "9 weeks")%>%mutate_if(is.character, as.factor)
  
  # 4th time point 
  
  t4status.sample <- t4[sample.int(nrow(t4), 500, replace=FALSE, prob=NULL),]
  t4status.sample <- t(t4status.sample)
  t4status.sample <- as.data.frame(cbind(rownames(t4status.sample), data.frame(t4status.sample, row.names=NULL)))
  t4status.sample <- t4status.sample %>% rename(CID=`rownames(t4status.sample)`)
  
  
  merge <- mergekkcca %>% select(dateN, CID, data.condition.column) %>%
    filter(dateN=="6 months")%>%
    right_join(t4status.sample, by="CID")
  
  merge <- merge[-which(is.na(merge$dateN)==T),]
  merge <- merge[!duplicated(merge[,"CID"]),]  
  
  # get the sensitivity/ specificity 
  for(i in 1:nrow(merge)){
    for(c in 4:ncol(merge)){
      if(merge[i,3]==1 & merge[i,c]==1){
        merge[i,c] <- "TP"
      } else if(merge[i,3]==1 & merge[i,c]==0){
        merge[i,c] <- "FP"
      } else if(merge[i,3]==0 & merge[i,c]==0){
        merge[i,c] <- "TN" 
      } else {
        merge[i,c] <- "FN" 
      }
    }
  }
  merge <- merge %>% mutate_if(is.character, as.factor)
  speclist <- list()
  
  # calculate how often each one occurs and get proportions 
  for(i in 4:ncol(merge)){
    speclist[[i-3]] <-  merge %>% group_by(merge[,i]) %>% tally() %>%
      mutate(proportion = n/sum(n))
  }
  
  test.specs<- list()
  for(i in 1:length(speclist)){
    test.specs[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
    colnames(test.specs[[i]]) <- c("TPR", "FPR")
    test.specs[[i]][1,1] <- speclist[[i]][4,2]/(speclist[[i]][4,2]+speclist[[i]][1,2])
    test.specs[[i]][1,2] <- speclist[[i]][2,2]/(speclist[[i]][3,2]+speclist[[i]][2,2])
  }
  
  roct4 <- rbindlist(test.specs)
  roct4 <- roct4 %>% mutate(time = "6 months")%>%mutate_if(is.character, as.factor)
  
  
  rocalltime <- bind_rows(roct1, roct2, roct3, roct4)
  return(rocalltime)
}




  