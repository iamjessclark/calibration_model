# Calibrating Schistosomiasis Diagnostics
# Jessica Clark, JM Prada, P.H.L.Lamberton
# Using the bayesian framework developed in the HMM to look at each time point independently
# Status will not be tracked across time points but will be estimated independently

# this is for just the KK data
CalModKK <- function (N, Ti, R, KK) {
  ## Set seed ##
  .RN1G.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  
  ## Initialize Status ##
  Status <- matrix(1,nrow=210, ncol=4)
  
  m <- "model {
    
    # Prior prevalence  #

    prev ~ dbeta(1,1)
    
    # Prior KKs #
    
    rtnb ~ dgamma(286.3235, 223.4650)
    
    for(t in 1:4){
      tkksh[t] ~ dgamma(83.88592,125.70331)
      tkkrt1[t] ~ dbeta(47.13542,633.08366)
      tkkrt[t] <- tkkrt1[t]/(1-tkkrt1[t])
    }
    
    # model 
    for(n in 1:N){	# run through pop
  
      for (t in 1:Ti){ # run through time
        Status[n,t] ~ dbern(prev)
        
        tKK[n,1,t] <- 0
        tKK[n,2,t] ~ dgamma(tkksh[t], tkkrt[t])
        
        lambda[n,t] <- tKK[n,Status[n,t]+1,t] 
          
          for(r in 1:R){ # run through repeated measures to set the baseline KK over the repeated measures, this is the LLH
             KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb)
          } # end or r loop
        } # end time
      } # end n loop
    
    #inits# .RNG.seed, .RNG.name, Status
    #data# N, Ti, R, KK
    #monitor# prev, rtnb, Status, tkksh, tkkrt
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
  prob <- array(rep(CCA,2),dim=c(N,Ti,2))
  
  m <- "model {
  # Prior prevalence #
  
  for(t in 1:4){
    prev[t] ~ dbeta(1,1)
  }
  
  # Prior KKs 
  rtnb ~ dgamma(286.3235, 223.4650)
  
  # make the mean sh/rt and sort them 
  # shapes
  
  #ssh[1] ~ dgamma(83.88592,125.70331)
  ssh[1] ~ dgamma(0.001,0.001)
  
  #ssh[2] ~ dgamma(20.57512,79.64904)
  ssh[2] ~ dgamma(0.001,0.001)
  
  ssh[3] ~ dgamma(0.001, 0.001)
  
  #ssh[4] ~ dgamma(53.59118,138.86460)
  ssh[4] ~ dgamma(0.001,0.001)
  
  # rates
  
  #srt1[1] ~ dbeta(47.13542,633.08366)
  srt1[1] ~ dbeta(1,1)
  srt[1] <- srt1[1]/(1-srt1[1])
  
  #srt2[2] ~ dbeta(6.946785,14.183304)
  srt2[2] ~ dbeta(1,1)
  srt[2] <- srt2[2]/(1-srt2[2])
  
  srt3[3] ~ dbeta(1,1)
  srt[3] <- srt3[3]/(1-srt3[3])
  
  #srt4[4] ~ dbeta(17.86595,77.49737)
  srt4[4] ~ dbeta(1,1)
  srt[4] <- srt4[4]/(1-srt4[4])
  
  ind <- order(ssh/srt)
  ptemp ~ dbern(0.5)
  ptemp1 ~ dbern(0.5)
  
  tkksh[1] <- ssh[ind[ptemp1+3]] # baseline 1 
  tkksh[2] <- ssh[ind[1+ptemp]] # 4w post - low
  tkksh[3] <- ssh[ind[1+(1-ptemp)]] # 4w post - low
  tkksh[4] <- ssh[ind[(1-ptemp1)+3]]
  
  tkkrt[1] <- srt[ind[ptemp1+3]] # baseline 1 
  tkkrt[2] <- srt[ind[1+ptemp]] # 4w post - low
  tkkrt[3] <- srt[ind[1+(1-ptemp)]] # 4w post - low
  tkkrt[4] <- srt[ind[(1-ptemp1)+3]] # baseline 2
  
  k ~ dgamma(0.001, 0.001)
  intercept ~ dnorm(0, 10^-6)
  
  # MODEL 
  for(n in 1:N){	# run through pop
    
    for (t in 1:Ti){ # run through timepoints
      
      Status[n,t] ~ dbern(prev[t])
      
      lambda[n,t] <- tKK[n,Status[n,t]+1,t] # this needs to go here because each time the r loop reiterates it will replace this value
        
        for( r in 1:R){  # run through repeat measures
          KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb) # generating the data with noise and then sampling from the dataset with the gamma sh/rt parameters?
      } # end or r loop

    
      tKK[n,1,t] <- 0                                                                                                                                                                                                                                                                                               
      tKK[n,2,t] ~ dgamma(tkksh[t], tkkrt[t])
      
      prob[n,t,1]~ dnorm(0,1.863282)T(0,4)
      prob[n,t,2] ~ dnorm(4 / (1 + exp(-k*(tKK[n,2,t]-intercept))), 3.47182)
      #3.093451  
      CCA[n,t] ~ dround(prob[n,t,Status[n,t]+1],0)
      
      } # end timestep t loop
    } #N
  
  #inits# .RNG.seed, .RNG.name, Status, prob   
  #data# N, Ti, R, KK, CCA
  #monitor#  rtnb, tkksh, tkkrt, prev, k, intercept, Status
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
  
  m <- "model {
  
  # MODEL 
   # Prior prevalence #
    
    for(t in 1:4){
      prev[t] ~ dbeta(1,1)
    }
   
    # Prior KKs 
  rtnb ~ dgamma(286.3235, 223.4650)
  
  # make the mean sh/rt and sort them 
  # shapes
  
  #ssh[1] ~ dgamma(83.88592,125.70331)
  ssh[1] ~ dgamma(0.001,0.001)
  
  #ssh[2] ~ dgamma(20.57512,79.64904)
  ssh[2] ~ dgamma(0.001,0.001)
  
  ssh[3] ~ dgamma(0.001, 0.001)
  
  #ssh[4] ~ dgamma(53.59118,138.86460)
  ssh[4] ~ dgamma(0.001,0.001)
  
  # rates
  
  #srt1[1] ~ dbeta(47.13542,633.08366)
  srt1[1] ~ dbeta(1,1)
  srt[1] <- srt1[1]/(1-srt1[1])
  
  #srt2[2] ~ dbeta(6.946785,14.183304)
  srt2[2] ~ dbeta(1,1)
  srt[2] <- srt2[2]/(1-srt2[2])
  
  srt3[3] ~ dbeta(1,1)
  srt[3] <- srt3[3]/(1-srt3[3])
  
  #srt4[4] ~ dbeta(17.86595,77.49737)
  srt4[4] ~ dbeta(1,1)
  srt[4] <- srt4[4]/(1-srt4[4])
  
  ind <- order(ssh/srt)
  ptemp ~ dbern(0.5)
  ptemp1 ~ dbern(0.5)
  
  tkksh[1] <- ssh[ind[ptemp1+3]] # baseline 1 
  tkksh[2] <- ssh[ind[1+ptemp]] # 4w post - low
  tkksh[3] <- ssh[ind[1+(1-ptemp)]] # 4w post - low
  tkksh[4] <- ssh[ind[(1-ptemp1)+3]]
  
  tkkrt[1] <- srt[ind[ptemp1+3]] # baseline 1 
  tkkrt[2] <- srt[ind[1+ptemp]] # 4w post - low
  tkkrt[3] <- srt[ind[1+(1-ptemp)]] # 4w post - low
  tkkrt[4] <- srt[ind[(1-ptemp1)+3]] # baseline 2
    
    k ~ dgamma(0.001, 0.001)
    intercept ~ dnorm(0, 10^-6)
    
    for(n in 1:N){	# run through pop
      
      for (t in 1:Ti){ # run through timepoints
      
            Status[n,t] ~ dbern(prev[t])
            
            lambda[n,t] <- tKK[n,Status[n,t]+1,t] 
            
            for( r in 1:R){  # run through repeat measures
              KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb) 
              } # end or r loop
      
            tKK[n,1,t] <- 0
            tKK[n,2,t] ~ dgamma(tkksh[t], tkkrt[t])
  
            prob[n,t,1]~ dnorm(0,1.245443)T(0,9)
            prob[n,t,2] ~ dnorm(9 / (1 + exp(-k*(tKK[n,2,t]-intercept))),1.245443)
            #1.093606
            
            CCA10[n,t] ~ dround(prob[n,t,Status[n,t]+1],0)
            
      } # t loop
    } # N loop   
    
 
  #inits# .RNG.seed, .RNG.name, Status, prob
  #data# N, Ti, R, KK, CCA10
  #monitor#  rtnb, tkksh, tkkrt, prev, k, intercept, Status
  
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
  return(unique(dt$child_id))
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
  return(unique(dt$CID))
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
  return(unique(dt$CID))
}

#### sort model output into time steps ####

time.steps <- function(model.output){
  
  t1 <- as.data.frame(model.output[,1:210])
  colnames(t1) <- CID
  t2 <- as.data.frame(model.output[,211:420])
  colnames(t2) <- CID
  t3 <- as.data.frame(model.output[,421:630])
  colnames(t3) <- CID
  t4 <- as.data.frame(model.output[,631:840])
  colnames(t4) <- CID
  
  t1[nrow(t1)+1,] <- colSums(t1)
  t1[nrow(t1),] <- t1[nrow(t1),]/40000
  
  t1.means <- as.data.frame(t1[nrow(t1),1:210])
  
  t1.means$time <-as.factor("pre-treatment")
  
  T1 <- reshape2::melt(t1.means)
  
  t2[nrow(t2)+1,] <- colSums(t2)
  t2[nrow(t2),] <- t2[nrow(t2),]/40000
  
  t2.means <- as.data.frame(t2[nrow(t2),1:210])
  
  t2.means$time <-as.factor("three-weeks")
  
  T2 <- reshape2::melt(t2.means)
  
  t3[nrow(t3)+1,] <- colSums(t3)
  t3[nrow(t3),] <- t3[nrow(t3),]/40000
  
  t3.means <- as.data.frame(t3[nrow(t3),1:210])
  
  t3.means$time <- as.factor("nine-weeks")
  
  T3 <- reshape2::melt(t3.means)
  
  t4[nrow(t4)+1,] <- colSums(t4)
  t4[nrow(t4),] <- t4[nrow(t4),]/40000
  
  t4.means <- as.data.frame(t4[nrow(t4),1:210])
  
  t4.means$time <- as.factor("six-months")
  
  T4 <- reshape2::melt(t4.means)
  
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

logcurvecca <- function(x,k,intercept){
  y <- 4 / (1 + exp(-k*(x-intercept)))
  return(y)
}


logcurvegscore <- function(x,k,intercept){
  y <- 9 / (1 + exp(-k*(x-intercept)))
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

#### ROC curve functions ####

## G-Score 

plotfuncgs <- function(fptpobj, aucobj){
  
  # auc values
  aucvals <- list()
  for(i in 1:500){
    aucvals[[i]] <- slot(aucobj, "y.values" )[[i]]
  }
  
  aucvals <- data.frame(do.call(rbind, aucvals))
  colnames(aucvals) <- "AUC"
  sdauc <- sd(aucvals$AUC)
  
  # TP/ FP values 
  tprate <- list()
  fprate <- list()
  for(i in 1:500){
    tprate[[i]] <- slot(fptpobj, "y.values")[[i]]
    fprate[[i]] <- slot(fptpobj, "x.values")[[i]]
  }
  
  tpraterate <- data.frame(do.call(rbind, tprate))
  tpraterate <- melt(tpraterate)
  colnames(tpraterate) <- c("cutoff", "tprate")
  fpratedf <- data.frame(do.call(rbind, fprate))
  fpratedf <- melt(fpratedf)
  colnames(fpratedf) <- c("cutoff", "fprate")
  fpratedf$cutoff <- dplyr::recode(fpratedf$cutoff, "X1"="Nothing is positive", 
                                   "X2"="G9","X3"="G8", "X4"="G7", 
                                   "X5"="G6", "X6"="G5", "X7"="G4",
                                   "X8"="G3", "X9"="G2", "X10"="G1")
  
  # make plot data
  plotdata <- data.frame(FP=fpratedf$fprate, TP=tpraterate$tprate, CUT=fpratedf$cutoff)
  
  
  plotdata <- melt(data=plotdata, id.vars = "CUT", measure.vars = c("TP", "FP"))
  plotdatasum <- plotdata %>%
    group_by(CUT, variable) %>%
    summarise(mean=mean(value), sd=sd(value))%>%
    pivot_wider(names_from = c(variable), values_from=c(mean, sd))%>%
    mutate(CUT=factor(CUT))
  
  # make plot 
  colours=c("#dadaeb","#bcbddc","#bfd3e6","#9ebcda","#74a9cf", "#8c96c6","#8c6bb1","#88419d","#810f7c", "black")
  plot <- ggplot()+
    geom_abline(intercept=0,slope=1) +
    geom_path(data=plotdatasum, aes(x=mean_FP, y=mean_TP), lwd=1) + 
    geom_point(data=plotdatasum, aes(x=mean_FP, y=mean_TP,colour=CUT), size=5, shape=1, stroke=3)+
    geom_errorbar(data=plotdatasum,aes(x=mean_FP, ymin=mean_TP-sd_TP, ymax=mean_TP+sd_TP))+
    geom_errorbarh(data=plotdatasum,aes(y=mean_TP, xmin=mean_FP-sd_FP, xmax=mean_FP+sd_FP))+
    annotate("text",x=0.97,y=0.02,label=paste("AUC = ",round(mean(aucvals$AUC),digits=2),"±", round(sdauc, 2),sep=""),hjust=1) +
    scale_colour_manual("Threshold Cutoff",values = rev(colours))+
    #scale_colour_gradientn("Threhsold Cutoff",colours=rainbow(14)[1:11], trans = 'reverse') +
    #annotate("text",x=0.97,y=0.06,label=paste("Sens = ",round(mean(specsen$sensitivity),digits=2),sep=""),hjust=1) +
    #annotate("text",x=0.97,y=0.08,label=paste("Spec = ",round(mean(specsen$specificity),digits=2),sep=""),hjust=1) 
    xlab("False Positive Rate") + ylab("True Positive Rate")+
    theme_bw()+
    theme(text = element_text(size=16))
  return(plot)
  
}

plotfuncgst4 <- function(fptpobj, aucobj){
  
  # auc values
  aucvals <- list()
  for(i in 1:500){
    aucvals[[i]] <- slot(aucobj, "y.values" )[[i]]
  }
  
  aucvals <- data.frame(do.call(rbind, aucvals))
  colnames(aucvals) <- "AUC"
  sdauc <- sd(aucvals$AUC)
  
  # TP/ FP values 
  tprate <- list()
  fprate <- list()
  for(i in 1:500){
    tprate[[i]] <- slot(fptpobj, "y.values")[[i]]
    fprate[[i]] <- slot(fptpobj, "x.values")[[i]]
  }
  
  tpraterate <- data.frame(do.call(rbind, tprate))
  tpraterate <- melt(tpraterate)
  colnames(tpraterate) <- c("cutoff", "tprate")
  fpratedf <- data.frame(do.call(rbind, fprate))
  fpratedf <- melt(fpratedf)
  colnames(fpratedf) <- c("cutoff", "fprate")
  fpratedf$cutoff <- dplyr::recode(fpratedf$cutoff, "X1"="Nothing is positive", "X2"="G10","X3"="G9", "X4"="G8", 
                                   "X5"="G7", "X6"="G6", "X7"="G5",
                                   "X8"="G4", "X9"="G3", "X10"="G2", "X11"="G1")

  # make plot data
  plotdata <- data.frame(FP=fpratedf$fprate, TP=tpraterate$tprate, CUT=fpratedf$cutoff)
  
  
  plotdata <- melt(data=plotdata, id.vars = "CUT", measure.vars = c("TP", "FP"))
  plotdatasum <- plotdata %>%
    group_by(CUT, variable) %>%
    summarise(mean=mean(value), sd=sd(value))%>%
    pivot_wider(names_from = c(variable), values_from=c(mean, sd))%>%
    mutate(CUT=factor(CUT))
  
  # make plot 
  colours=c("#dadaeb","#bcbddc","#bfd3e6","#9ebcda","#74a9cf", "#8c96c6","#8c6bb1","#88419d","#810f7c","#4d004b", "black")
  plot <- ggplot()+
    geom_abline(intercept=0,slope=1) +
    geom_path(data=plotdatasum, aes(x=mean_FP, y=mean_TP), lwd=1) + 
    geom_point(data=plotdatasum, aes(x=mean_FP, y=mean_TP,colour=CUT), size=5, shape=1, stroke=3)+
    geom_errorbar(data=plotdatasum,aes(x=mean_FP, ymin=mean_TP-sd_TP, ymax=mean_TP+sd_TP))+
    geom_errorbarh(data=plotdatasum,aes(y=mean_TP, xmin=mean_FP-sd_FP, xmax=mean_FP+sd_FP))+
    annotate("text",x=0.97,y=0.02,label=paste("AUC = ",round(mean(aucvals$AUC),digits=2),"±", round(sdauc, 2),sep=""),hjust=1) +
    scale_colour_manual("Threshold Cutoff",values = rev(colours))+
    #scale_colour_gradientn("Threhsold Cutoff",colours=rainbow(14)[1:11], trans = 'reverse') +
    #annotate("text",x=0.97,y=0.06,label=paste("Sens = ",round(mean(specsen$sensitivity),digits=2),sep=""),hjust=1) +
    #annotate("text",x=0.97,y=0.08,label=paste("Spec = ",round(mean(specsen$specificity),digits=2),sep=""),hjust=1) 
    xlab("False Positive Rate") + ylab("True Positive Rate")+
    theme_bw()+
    theme(text = element_text(size=16))
  return(plot)
  
}


## CCA 
plotfunccca <- function(fptpobj, aucobj){
  
  # auc values
  aucvals <- list()
  for(i in 1:500){
    aucvals[[i]] <- slot(aucobj, "y.values" )[[i]]
  }
  
  aucvals <- data.frame(do.call(rbind, aucvals))
  colnames(aucvals) <- "AUC"
  sdauc <- sd(aucvals$AUC)
  
  # TP/ FP values 
  tprate <- list()
  fprate <- list()
  for(i in 1:500){
    tprate[[i]] <- slot(fptpobj, "y.values")[[i]]
    fprate[[i]] <- slot(fptpobj, "x.values")[[i]]
  }
  
  tpraterate <- data.frame(do.call(rbind, tprate))
  tpraterate <- melt(tpraterate)
  colnames(tpraterate) <- c("cutoff", "tprate")
  fpratedf <- data.frame(do.call(rbind, fprate))
  fpratedf <- melt(fpratedf)
  colnames(fpratedf) <- c("cutoff", "fprate")
  fpratedf$cutoff <- dplyr::recode(fpratedf$cutoff, "X1"="Nothing is positive", "X2"="+++","X3"="++", "X4"="+", "X5"="Trace", "X6"="Negative")
  #fpratedf$cutoff <- as.numeric(as.character(fpratedf$cutoff))
  
  # make plot data
  plotdata <- data.frame(FP=fpratedf$fprate, TP=tpraterate$tprate, CUT=fpratedf$cutoff)
  
  
  plotdata <- melt(data=plotdata, id.vars = "CUT", measure.vars = c("TP", "FP"))
  plotdatasum <- plotdata %>%
    group_by(CUT, variable) %>%
    summarise(mean=mean(value), sd=sd(value))%>%
    pivot_wider(names_from = c(variable), values_from=c(mean, sd))%>%
    mutate(CUT=factor(CUT))
  
  # make plot 
  colours=c( "#fdd0a2","#fd8d3c","#f16913","#d94801","#8c2d04",  "black")
  plot <- ggplot()+
    geom_abline(intercept=0,slope=1) +
    geom_path(data=plotdatasum, aes(x=mean_FP, y=mean_TP), lwd=1) + 
    geom_point(data=plotdatasum, aes(x=mean_FP, y=mean_TP,colour=CUT), size=5, shape=1, stroke=3)+
    geom_errorbar(data=plotdatasum,aes(x=mean_FP, ymin=mean_TP-sd_TP, ymax=mean_TP+sd_TP))+
    geom_errorbarh(data=plotdatasum,aes(y=mean_TP, xmin=mean_FP-sd_FP, xmax=mean_FP+sd_FP))+
    annotate("text",x=0.97,y=0.02,label=paste("AUC = ",round(mean(aucvals$AUC),digits=2),"±", round(sdauc, 2),sep=""),hjust=1) +
    scale_colour_manual("Threshold Cutoff",values = rev(colours))+
    #scale_colour_gradientn("Threhsold Cutoff",colours=rainbow(14)[1:11], trans = 'reverse') +
    #annotate("text",x=0.97,y=0.06,label=paste("Sens = ",round(mean(specsen$sensitivity),digits=2),sep=""),hjust=1) +
    #annotate("text",x=0.97,y=0.08,label=paste("Spec = ",round(mean(specsen$specificity),digits=2),sep=""),hjust=1) 
    xlab("False Positive Rate") + ylab("True Positive Rate")+
    theme_bw()+
    theme(text = element_text(size=16))
  return(plot)
  
}



#### roc curve function ####

diag_specgscoret1 <- function(t1, condition.data.frame){

  # join to the condition
  merge <- condition.data.frame %>%
    right_join(t1, by="CID")
  
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
  
  # calculate how often each one occurs and get proportions 
  speclist <- list()
  
  for(i in 4:ncol(merge)){
    merge[,i] <- factor(merge[,i], levels=c("TP", "FP", "FN", "TN"))
  }
  
  for(i in 4:ncol(merge)){
    speclist[[i-3]] <- table(merge[,i])
  }
  
  test.specs<- list()
  for(i in 1:length(speclist)){
    test.specs[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
    colnames(test.specs[[i]]) <- c("TPR", "FPR")
    
    test.specs[[i]][1,1] <- speclist[[i]][[1]]/(speclist[[i]][[1]]+speclist[[i]][[3]])
    test.specs[[i]][1,2] <- speclist[[i]][[2]]/(speclist[[i]][[2]]+speclist[[i]][[4]])
  }
  
  roct1 <- rbindlist(test.specs)
  roct1 <- roct1 %>% mutate(time = "pre-T")%>% mutate_if(is.character, as.factor)

   return(roct1)
}

diag_specgscoret2 <- function(t2, condition.data.frame ){
  # 2nd time point 
  
  merge <- condition.data.frame %>%
    right_join(t2, by="CID")%>%
    filter(dateN=="3 weeks")
  
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
  
  # calculate how often each one occurs and get proportions 
  speclist <- list()
  
  for(i in 4:ncol(merge)){
    merge[,i] <- factor(merge[,i], levels=c("TP", "FP", "FN", "TN"))
  }
  
  for(i in 4:ncol(merge)){
    speclist[[i-3]] <- table(merge[,i])
  }
  
  test.specs<- list()
  for(i in 1:length(speclist)){
    test.specs[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
    colnames(test.specs[[i]]) <- c("TPR", "FPR")
    
    test.specs[[i]][1,1] <- speclist[[i]][[1]]/(speclist[[i]][[1]]+speclist[[i]][[3]])
    test.specs[[i]][1,2] <- speclist[[i]][[2]]/(speclist[[i]][[2]]+speclist[[i]][[4]])
  }

  roct2 <- rbindlist(test.specs)
  roct2 <- roct2 %>% mutate(time = "3 weeks")%>%mutate_if(is.character, as.factor)
  
  return(roct2)
  
}

diag_specgscoret3 <- function(t3, condition.data.frame ){
  
  merge <- condition.data.frame %>%
    right_join(t3, by="CID")%>%
    filter(dateN=="9 weeks")
  
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
  
  speclist <- list()
  
  for(i in 4:ncol(merge)){
    merge[,i] <- factor(merge[,i], levels=c("TP", "FP", "FN", "TN"))
  }
  
  for(i in 4:ncol(merge)){
    speclist[[i-3]] <- table(merge[,i])
  }
  
  test.specs<- list()
  for(i in 1:length(speclist)){
    test.specs[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
    colnames(test.specs[[i]]) <- c("TPR", "FPR")
    
    test.specs[[i]][1,1] <- speclist[[i]][[1]]/(speclist[[i]][[1]]+speclist[[i]][[3]])
    test.specs[[i]][1,2] <- speclist[[i]][[2]]/(speclist[[i]][[2]]+speclist[[i]][[4]])
  }
  
  roct3 <- rbindlist(test.specs)
  roct3 <- roct3 %>% mutate(time = "9 weeks")%>%mutate_if(is.character, as.factor)
  
  return(roct3)
  
}
  # 4th time point 
diag_specgscoret4 <- function(t4, condition.data.frame){
  
  merge <- condition.data.frame %>%
    right_join(t4, by="CID")%>%
    filter(dateN=="6 months")
  
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
  
  # get the sensitivity/ specificity 
  speclist <- list()
  
  for(i in 4:ncol(merge)){
    merge[,i] <- factor(merge[,i], levels=c("TP", "FP", "FN", "TN"))
  }
  
  for(i in 4:ncol(merge)){
    speclist[[i-3]] <- table(merge[,i])
  }
  
  test.specs<- list()
  for(i in 1:length(speclist)){
    test.specs[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
    colnames(test.specs[[i]]) <- c("TPR", "FPR")
    
    test.specs[[i]][1,1] <- speclist[[i]][[1]]/(speclist[[i]][[1]]+speclist[[i]][[3]])
    test.specs[[i]][1,2] <- speclist[[i]][[2]]/(speclist[[i]][[2]]+speclist[[i]][[4]])
  }
  
  roct4 <- rbindlist(test.specs)
  roct4 <- roct4 %>% mutate(time = "6 months")%>%mutate_if(is.character, as.factor)
  
  return(roct4)
}

# calculating the optimal point on the ROC curve by maximising Youden's index
# this is when you maximise the sum of the spec and sens
# according to Kaivanto, 2008, Journal of Clinical Epi this is fine if the aim is maximise correct ID of disease 
# and if the cost of misclassifying the diseased is equal to mis-classifying the non-diseased. 
# I think this is justiable. 
rocr_sensspec <- function(x, class) {
  pred <- ROCR::prediction(x, class)
  perf <- ROCR::performance(pred, "sens", "spec")
  sens <- slot(perf, "y.values")[[1]]
  spec <- slot(perf, "x.values")[[1]]
  cut <- slot(perf, "alpha.values")[[1]]
  cut[which.max(sens + spec)]
}

# extracting the kk data from the array 
getKKtime <- function(ArrayName, TimeStep){
  kkTime <- as.data.frame(cbind(ArrayName[,,1][,TimeStep],ArrayName[,,2][,TimeStep],ArrayName[,,3][,TimeStep],ArrayName[,,4][,TimeStep],ArrayName[,,5][,TimeStep],ArrayName[,,6][,TimeStep]))%>%
    mutate(CID=CID, 
           meancount=rowMeans(., na.rm = T),
           kkstatus=case_when(meancount>0~1,
                              meancount==0~0, 
                              TRUE~NA_real_))
  return(kkTime)
}

  