require(data.table)
require(reshape2)
require(tidyverse)
require(runjags)
require(rjags)
require(coda)
require(viridis)
require(ggpubr)
require(ROCR)

require(truncnorm)

set.seed(2)
# load data #
source("diagnositic.calibration.model.R")

kkdata <- read.csv("KK.4tp.csv")
ccadata <- read.csv("cca.clean.csv")

kk <- loadKKdata("KK.4tp.csv")
KKIDs <- getKKChildIDs("KK.4tp.csv")

cca <- loadCCAdata("cca.clean.csv")
CCAIDs <- getCCAChildIDs("cca.clean.csv")
CID <- as.character(CCAIDs[(CCAIDs %in% KKIDs)])

kk <- kk[match(CID, KKIDs),,]
cca <- cca[match(CID, CCAIDs),]

ccagscore <- loadgscoredata("cca.clean.csv")
gscoreCCAIDS <- getCCAChildIDsgscore("cca.clean.csv") 
gscoreCID <- as.character(gscoreCCAIDS[(gscoreCCAIDS %in% KKIDs)])
ccagscore <- ccagscore[match(gscoreCID,gscoreCCAIDS),]

#### Getting scores and kk counts from the data that goes into the model ####
# kato katz

getKKtime <- function(ArrayName, TimeStep){
  kkTime <- as.data.frame(cbind(ArrayName[,,1][,TimeStep],ArrayName[,,2][,TimeStep],ArrayName[,,3][,TimeStep],ArrayName[,,4][,TimeStep],ArrayName[,,5][,TimeStep],ArrayName[,,6][,TimeStep]))%>%
    mutate(CID=CID, 
           meancount=rowMeans(., na.rm = T),
           kkstatus=case_when(meancount>0~1,
                              meancount==0~0, 
                              TRUE~NA_real_))
  return(kkTime)
}

kkpre <- as.data.frame(cbind(kk[,,1][,1],kk[,,2][,1],kk[,,3][,1],kk[,,4][,1],kk[,,5][,1],kk[,,6][,1]))%>%
  mutate(CID=CID, 
         meancount=rowMeans(., na.rm = T),
         kkstatus=case_when(meancount>0~1,
                            meancount==0~0, 
                            TRUE~NA_real_))


kkthreewk <- as.data.frame(cbind(kk[,,1][,2],kk[,,2][,2],kk[,,3][,2],kk[,,4][,2],kk[,,5][,2],kk[,,6][,2]))%>%
  mutate(CID=CID, 
         meancount=rowMeans(., na.rm = T),
         kkstatus=case_when(meancount>0~1,
                            meancount==0~0, 
                            TRUE~NA_real_))

kkninewk <- as.data.frame(cbind(kk[,,1][,3],kk[,,2][,3],kk[,,3][,3],kk[,,4][,3],kk[,,5][,3],kk[,,6][,3]))%>%
  mutate(CID=CID, 
         meancount=rowMeans(., na.rm = T),
         kkstatus=case_when(meancount>0~1,
                            meancount==0~0, 
                            TRUE~NA_real_))

kksixmth <- as.data.frame(cbind(kk[,,1][,4],kk[,,2][,4],kk[,,3][,4],kk[,,4][,4],kk[,,5][,4],kk[,,6][,4] ))%>%
  mutate(CID=CID, 
         meancount=rowMeans(., na.rm = T),
         kkstatus=case_when(meancount>0~1,
                            meancount==0~0, 
                            TRUE~NA_real_))
# cca data 

ccascores <- as.data.frame(cca)
colnames(ccascores) <- c("Pre-T", "3 weeks", "9 weeks", "6 months")
ccascores$CID <- CID
ccascores <- melt(ccascores, id.vars = c("CID"), measure.vars = c("Pre-T", "3 weeks", "9 weeks", "6 months"))
colnames(ccascores) <- c("CID", "time", "CCA")

ccapre <- subset(ccascores, time=="Pre-T")
ccathreewk <- subset(ccascores, time=="3 weeks")
ccaninewk <- subset(ccascores, time=="9 weeks")
ccasixmth <- subset(ccascores, time=="6 months")

# gscore 

gscoredata <- as.data.frame(ccagscore)
colnames(gscoredata) <- c("Pre-T", "3 weeks", "9 weeks", "6 months")
gscoredata$CID <- CID
gscoredata <- melt(gscoredata, id.vars = c("CID"), measure.vars = c("Pre-T", "3 weeks", "9 weeks", "6 months"))
colnames(gscoredata) <- c("CID", "time", "Gscore")

gscorepre <- subset(gscoredata, time=="Pre-T")
gscorethreewk <- subset(gscoredata, time=="3 weeks")
gscoreninewk <- subset(gscoredata, time=="9 weeks")
gscoresixmth <- subset(gscoredata, time=="6 months")


predata <- merge(ccapre, kkpre, by="CID")
predata <- merge(gscorepre, predata, by="CID")%>%
  select(CID, Gscore, CCA, kkstatus)%>%
  mutate(ccanegpos = 1,
         ccaTpos = ifelse(CCA == 1 |CCA ==  2 | CCA == 3 |CCA == 4, 1,0),
         cca1pos = ifelse(CCA ==  2 | CCA == 3 |CCA == 4, 1,0),
         cca2pos = ifelse(CCA == 3 |CCA == 4, 1,0),
         cca3pos = ifelse(CCA == 4, 1,0),
         gscore1pos = ifelse(Gscore>=1, 1, 0), 
         gscore2pos = ifelse(Gscore>=2, 1, 0), 
         gscore3pos = ifelse(Gscore>=3, 1, 0), 
         gscore4pos = ifelse(Gscore>=4, 1, 0), 
         gscore5pos = ifelse(Gscore>=5, 1, 0), 
         gscore6pos = ifelse(Gscore>=6, 1, 0), 
         gscore7pos = ifelse(Gscore>=7, 1, 0), 
         gscore8pos = ifelse(Gscore>=8, 1, 0),
         gscore9pos = ifelse(Gscore>=9, 1, 0), 
         gscore10pos = ifelse(Gscore==10, 1, 0))

threewkdata <- merge(ccathreewk, kkthreewk, by="CID")
threewkdata <- merge(gscorethreewk, threewkdata, by="CID")%>%
  select(CID, Gscore, CCA, kkstatus)%>%
  mutate(ccanegpos = 1,
         ccaTpos = ifelse(CCA == 1 |CCA ==  2 | CCA == 3 |CCA == 4, 1,0),
         cca1pos = ifelse(CCA ==  2 | CCA == 3 |CCA == 4, 1,0),
         cca2pos = ifelse(CCA == 3 |CCA == 4, 1,0),
         cca3pos = ifelse(CCA == 4, 1,0),
         gscore1pos = ifelse(Gscore>=1, 1, 0), 
         gscore2pos = ifelse(Gscore>=2, 1, 0), 
         gscore3pos = ifelse(Gscore>=3, 1, 0), 
         gscore4pos = ifelse(Gscore>=4, 1, 0), 
         gscore5pos = ifelse(Gscore>=5, 1, 0), 
         gscore6pos = ifelse(Gscore>=6, 1, 0), 
         gscore7pos = ifelse(Gscore>=7, 1, 0), 
         gscore8pos = ifelse(Gscore>=8, 1, 0),
         gscore9pos = ifelse(Gscore>=9, 1, 0), 
         gscore10pos = ifelse(Gscore==10, 1, 0))


ninewkdata <- merge(ccaninewk, kkninewk, by="CID")
ninewkdata <- merge(gscoreninewk, ninewkdata, by="CID")%>%
  select(CID, Gscore, CCA, kkstatus)%>%
  mutate(ccanegpos = 1,
         ccaTpos = ifelse(CCA == 1 |CCA ==  2 | CCA == 3 |CCA == 4, 1,0),
         cca1pos = ifelse(CCA ==  2 | CCA == 3 |CCA == 4, 1,0),
         cca2pos = ifelse(CCA == 3 |CCA == 4, 1,0),
         cca3pos = ifelse(CCA == 4, 1,0),
         gscore1pos = ifelse(Gscore>=1, 1, 0), 
         gscore2pos = ifelse(Gscore>=2, 1, 0), 
         gscore3pos = ifelse(Gscore>=3, 1, 0), 
         gscore4pos = ifelse(Gscore>=4, 1, 0), 
         gscore5pos = ifelse(Gscore>=5, 1, 0), 
         gscore6pos = ifelse(Gscore>=6, 1, 0), 
         gscore7pos = ifelse(Gscore>=7, 1, 0), 
         gscore8pos = ifelse(Gscore>=8, 1, 0),
         gscore9pos = ifelse(Gscore>=9, 1, 0), 
         gscore10pos = ifelse(Gscore==10, 1, 0))

sixmonthdata <- merge(ccasixmth, kksixmth, by="CID")
sixmonthdata <- merge(gscoresixmth, sixmonthdata, by="CID")%>%
  select(CID, Gscore, CCA, kkstatus)%>%
  mutate(ccanegpos = 1,
         ccaTpos = ifelse(CCA == 1 |CCA ==  2 | CCA == 3 |CCA == 4, 1,0),
         cca1pos = ifelse(CCA ==  2 | CCA == 3 |CCA == 4, 1,0),
         cca2pos = ifelse(CCA == 3 |CCA == 4, 1,0),
         cca3pos = ifelse(CCA == 4, 1,0),
         gscore1pos = ifelse(Gscore>=1, 1, 0), 
         gscore2pos = ifelse(Gscore>=2, 1, 0), 
         gscore3pos = ifelse(Gscore>=3, 1, 0), 
         gscore4pos = ifelse(Gscore>=4, 1, 0), 
         gscore5pos = ifelse(Gscore>=5, 1, 0), 
         gscore6pos = ifelse(Gscore>=6, 1, 0), 
         gscore7pos = ifelse(Gscore>=7, 1, 0), 
         gscore8pos = ifelse(Gscore>=8, 1, 0),
         gscore9pos = ifelse(Gscore>=9, 1, 0), 
         gscore10pos = ifelse(Gscore==10, 1, 0))

#### KK and CCA Model ####

kkcca.model <- CalModKKCCA(210,4,6,kk,cca)
prevKKcca <- density(c(kkcca.model$mcmc[[1]][,"prev[1]"], kkcca.model$mcmc[[2]][,"prev[1]"]))
prevKKcca <- density(c(kkcca.model$mcmc[[1]][,"prev[2]"], kkcca.model$mcmc[[2]][,"prev[2]"]))
prevKKcca <- density(c(kkcca.model$mcmc[[1]][,"prev[3]"], kkcca.model$mcmc[[2]][,"prev[3]"]))
prevKKcca <- density(c(kkcca.model$mcmc[[1]][,"prev[4]"], kkcca.model$mcmc[[2]][,"prev[4]"]))

rtnbKKcca <- density(c(kkcca.model$mcmc[[1]][,"rtnb"], kkcca.model$mcmc[[2]][,"rtnb"]))
klogKKcca <- density(c(kkcca.model$mcmc[[1]][,"k"], kkcca.model$mcmc[[2]][,"k"]))
interceptKKcca <- density(c(kkcca.model$mcmc[[1]][,"intercept"], kkcca.model$mcmc[[2]][,"intercept"]))

sht1KKcca <- density(c(kkcca.model$mcmc[[1]][,"tkksh[1]"], kkcca.model$mcmc[[2]][,"tkksh[1]"]))
sht2KKcca <- density(c(kkcca.model$mcmc[[1]][,"tkksh[2]"], kkcca.model$mcmc[[2]][,"tkksh[2]"]))
sht3KKcca <- density(c(kkcca.model$mcmc[[1]][,"tkksh[3]"], kkcca.model$mcmc[[2]][,"tkksh[3]"]))
sht4KKcca <- density(c(kkcca.model$mcmc[[1]][,"tkksh[4]"], kkcca.model$mcmc[[2]][,"tkksh[4]"]))

rtt1KKcca <- density(c(kkcca.model$mcmc[[1]][,"tkkrt[1]"], kkcca.model$mcmc[[2]][,"tkkrt[1]"]))
rtt2KKcca <- density(c(kkcca.model$mcmc[[1]][,"tkkrt[2]"], kkcca.model$mcmc[[2]][,"tkkrt[2]"]))
rtt3KKcca <- density(c(kkcca.model$mcmc[[1]][,"tkkrt[3]"], kkcca.model$mcmc[[2]][,"tkkrt[3]"]))
rtt4KKcca <- density(c(kkcca.model$mcmc[[1]][,"tkkrt[4]"], kkcca.model$mcmc[[2]][,"tkkrt[4]"]))

model.outputkkcca <- as.data.frame(as.matrix(as.mcmc(kkcca.model)))
status.kkcca <- model.outputkkcca[,16:ncol(model.outputkkcca)]
# sample for roc curve
cca.status.sample <- status.kkcca[sample.int(nrow(status.kkcca), 500, replace=FALSE, prob=NULL),]

# random sampling for the log curve
kintkkcca <- as.data.frame(cbind(c(kkcca.model$mcmc[[1]][,"k"], kkcca.model$mcmc[[2]][,"k"]),c(kkcca.model$mcmc[[1]][,"intercept"], kkcca.model$mcmc[[2]][,"intercept"])))
kkcca.sample <- kintkkcca[sample.int(nrow(kintkkcca), 500, replace=FALSE, prob=NULL),]

# make prob infcted by time status object 
status.time.KKcca <- time.steps(status.kkcca)

anno <- data.frame(c("Baseline", "ThreeWeeks", "NineWeeks", "SixMonths"))
colnames(anno) <- "time"
anno$xcoord <- 0
anno$ycoord <-200
anno$label <- as.factor(c("A", "B", "C", "D"))


# probability of infection at each time point #
pdf("probinfKKcca.pdf",width=3,height=8)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfrow=c(1,1))

prob.inf.kkcca <- qplot(value, data = status.time.KKcca, geom = "histogram", fill=time, binwidth=0.02) +
  scale_fill_manual(values = c("#e7298a", "#1b9e77", "#d95f02", "#7570b3"))+
  geom_text(data=anno, aes(x = xcoord, y = ycoord, label = label))+
  facet_grid(time~.)+
  theme_bw()+coord_cartesian(ylim=c(0, 200), xlim=c(0,1))+theme(legend.position = "none")+
  ylab("Number of Hosts") + xlab("Probability of Being Infected")

prob.inf.kkcca
dev.off()

t1cca <- as.data.frame(cca.status.sample [,1:210])
t2cca <- as.data.frame(cca.status.sample[,211:420])
t3cca <- as.data.frame(cca.status.sample[,421:630])
t4cca <- as.data.frame(cca.status.sample[,631:840])

colnames(t1cca) <- CID
colnames(t2cca) <- CID
colnames(t3cca) <- CID
colnames(t4cca) <- CID

t1cca <- t(t1cca)
t2cca <- t(t2cca)
t3cca <- t(t3cca)
t4cca <- t(t4cca)

t1cca <- as.data.frame(cbind(rownames(t1cca), data.frame(t1cca, row.names=NULL)))
t1cca <- t1cca %>% rename(CID=`rownames(t1cca)`)

t2cca <- as.data.frame(cbind(rownames(t2cca), data.frame(t2cca, row.names=NULL)))
t2cca <- t2cca %>% rename(CID=`rownames(t2cca)`)

t3cca <- as.data.frame(cbind(rownames(t3cca), data.frame(t3cca, row.names=NULL)))
t3cca <- t3cca %>% rename(CID=`rownames(t3cca)`)

t4cca <- as.data.frame(cbind(rownames(t4cca), data.frame(t4cca, row.names=NULL)))
t4cca <- t4cca %>% rename(CID=`rownames(t4cca)`)

#### CCA ROCR ####

# pre treatment 

predatacca <- predata %>% dplyr::select(CID, CCA) %>%
  right_join(t1, by="CID")
predatacca <- predatacca[-which(is.na(predatacca$CCA)==T),]
predictionspreT <- rep(list(predatacca$CCA), 500)
labelslistpreT <- list()
for(i in 3:ncol(predatacca)){
  labelslistpreT[[i-2]] <- predatacca[,i]
}

# three weeks 

threeweekscca <- threewkdata %>% dplyr::select(CID, CCA) %>%
  right_join(t2, by="CID")
threeweekscca <- threeweekscca[-which(is.na(threeweekscca$CCA)==T),]
predictionsthreewk <- rep(list(threeweekscca$CCA), 500)
labelslistthreewk <- list()
for(i in 3:ncol(threeweekscca)){
  labelslistthreewk[[i-2]] <- threeweekscca[,i]
}

# nine weeks 

nineweekscca <- ninewkdata %>% dplyr::select(CID, CCA) %>%
  right_join(t3, by="CID")
nineweekscca <- nineweekscca[-which(is.na(nineweekscca$CCA)==T),]
predictionsninewk <- rep(list(nineweekscca$CCA), 500)
labelslistninewk <- list()
for(i in 3:ncol(nineweekscca)){
  labelslistninewk[[i-2]] <- nineweekscca[,i]
}


# six months 

sixmonthscca <- sixmonthdata %>% dplyr::select(CID, CCA) %>%
  right_join(t4, by="CID")
sixmonthscca <- sixmonthscca[-which(is.na(sixmonthscca$CCA)==T),]
predictionssixmth <- rep(list(sixmonthscca$CCA), 500)
labelslistsixmonth <- list()
for(i in 3:ncol(sixmonthscca)){
  labelslistsixmonth[[i-2]] <- sixmonthscca[,i]
}

# plot 
sort(table(cutoff2),decreasing=TRUE)[1:2]
# baseline
manypredt1 <- prediction(predictionspreT, labelslistpreT)
auc.perfccat1 <- performance(manypredt1, measure = "auc")
sensspec.perfccat1 <- performance(manypredt1, "sens", "spec")
auccca1 <- t(rbindlist(list(auc.perfccat1@y.values)))
quantile(auccca1,  probs = c(0.025,0.5, 0.975))

sens <- list()
spec <- list()
cut <- list()
cutoff <- vector()
for(i in 1:500){
  sens[[i]] <- slot(sensspec.perfccat1, "y.values")[[i]]
  spec[[i]] <- slot(sensspec.perfccat1, "x.values")[[i]]
  cut[[i]] <- slot(sensspec.perfccat1, "alpha.values")[[i]]
  
  cutoff[i] <- cut[[i]][which.max(sens[[i]] + spec[[i]])]
}


many.roc.ccat1 <-  performance(manypredt1, measure = "tpr", x.measure = "fpr")
#pdf("pretTroccca.pdf") 
plot(many.roc.ccat1, col="grey82", lty=3)
plot(many.roc.ccat1, lwd=3,average="threshold", add=TRUE)
abline(a=0, b= 1)
dev.off()

# three weeks 
manypredt2 <- prediction(predictionsthreewk, labelslistthreewk)
auc.perfccat2 <-  performance(manypredt2, measure = "auc")
sensspec.perfccat2 <- performance(manypredt2, "sens", "spec")
auccca2 <- t(rbindlist(list(auc.perfccat2@y.values)))
quantile(auccca2,  probs = c(0.025,0.5, 0.975))

sens <- list()
spec <- list()
cut <- list()
cutoff2 <- vector()
for(i in 1:500){
  sens[[i]] <- slot(sensspec.perfccat2, "y.values")[[i]]
  spec[[i]] <- slot(sensspec.perfccat2, "x.values")[[i]]
  cut[[i]] <- slot(sensspec.perfccat2, "alpha.values")[[i]]
  
  cutoff2[i] <- cut[[i]][which.max(sens[[i]] + spec[[i]])]
}

many.roc.ccat2 <-  performance(manypredt2, measure = "tpr", x.measure = "fpr")
#pdf("threewkroccca.pdf") 
plot(many.roc.ccat2, col="grey82", lty=3)
plot(many.roc.ccat2, lwd=3,  add=TRUE)
abline(a=0, b= 1)
dev.off()

# nine weeks 
manypredt3 <- prediction(predictionsninewk, labelslistninewk)
auc.perfccat3 <-  performance(manypredt3, measure = "auc")
sensspec.perfccat3 <- performance(manypredt3, "sens", "spec")
auccca3 <- t(rbindlist(list(auc.perfccat3@y.values)))
quantile(auccca3,  probs = c(0.025,0.5, 0.975))

sens <- list()
spec <- list()
cut <- list()
cutoff3 <- vector()
for(i in 1:500){
  sens[[i]] <- slot(sensspec.perfccat3, "y.values")[[i]]
  spec[[i]] <- slot(sensspec.perfccat3, "x.values")[[i]]
  cut[[i]] <- slot(sensspec.perfccat3, "alpha.values")[[i]]
  
  cutoff3[i] <- cut[[i]][which.max(sens[[i]] + spec[[i]])]
}

many.roc.ccat3 <-  performance(manypredt3, measure = "tpr", x.measure = "fpr")
pdf("ninewkroccca.pdf") 
plot(many.roc.ccat3, col="grey82", lty=3)
plot(many.roc.ccat3, lwd=3, avg="vertical", spread.estimate="boxplot", add=TRUE)
abline(a=0, b= 1)
dev.off()

manypredt4 <- prediction(predictionssixmth, labelslistsixmonth)
auc.perfccat4 <-  performance(manypredt4, measure = "auc")
sensspec.perfccat4 <- performance(manypredt4, "sens", "spec")
auccca4 <- t(rbindlist(list(auc.perfccat4@y.values)))
quantile(auccca4,  probs = c(0.025,0.5, 0.975))

sens <- list()
spec <- list()
cut <- list()
cutoff4 <- vector()
for(i in 1:500){
  sens[[i]] <- slot(sensspec.perfccat4, "y.values")[[i]]
  spec[[i]] <- slot(sensspec.perfccat4, "x.values")[[i]]
  cut[[i]] <- slot(sensspec.perfccat4, "alpha.values")[[i]]
  
  cutoff4[i] <- cut[[i]][which.max(sens[[i]] + spec[[i]])]
}


many.roc.ccat4 <-  performance(manypredt4, measure = "tpr", x.measure = "fpr")
pdf("sixmthroccca.pdf") 
plot(many.roc.ccat4, col="grey82", lty=3)
plot(many.roc.ccat4, lwd=3, avg="vertical", spread.estimate="boxplot", add=TRUE)
abline(a=0, b= 1)
dev.off()

#### KK and Gscore Model ####
kkGscore.model <- CalModKKGScore(210,4,6,kk,ccagscore)

prevKKGScore <- density(c(kkGscore.model$mcmc[[1]][,"prev"], kkGscore.model$mcmc[[2]][,"prev"]))

rtnbKKgscore <- density(c(kkGscore.model$mcmc[[1]][,"rtnb"], kkGscore.model$mcmc[[2]][,"rtnb"]))

sht1KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkksh[1]"], kkGscore.model$mcmc[[2]][,"tkksh[1]"]))
sht2KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkksh[2]"], kkGscore.model$mcmc[[2]][,"tkksh[2]"]))
sht3KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkksh[3]"], kkGscore.model$mcmc[[2]][,"tkksh[3]"]))
sht4KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkksh[4]"], kkGscore.model$mcmc[[2]][,"tkksh[4]"]))

rtt1KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkkrt[1]"], kkGscore.model$mcmc[[2]][,"tkkrt[1]"]))
rtt2KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkkrt[2]"], kkGscore.model$mcmc[[2]][,"tkkrt[2]"]))
rtt3KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkkrt[3]"], kkGscore.model$mcmc[[2]][,"tkkrt[3]"]))
rtt4KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkkrt[4]"], kkGscore.model$mcmc[[2]][,"tkkrt[4]"]))

klogKKgscore <- (c(kkGscore.model$mcmc[[1]][,"k"], kkGscore.model$mcmc[[2]][,"k"]))
interceptKKgscore <- (c(kkGscore.model$mcmc[[1]][,"intercept"], kkGscore.model$mcmc[[2]][,"intercept"]))

outputkkGscore <- as.data.frame(as.matrix(as.mcmc(kkGscore.model)))
status.kkgscore <- outputkkGscore[,13:ncol(outputkkGscore)]

#random samples to make CI from

kintkkgscore <- as.data.frame(cbind(c(kkGscore.model$mcmc[[1]][,"k"], kkGscore.model$mcmc[[2]][,"k"]),c(kkGscore.model$mcmc[[1]][,"intercept"], kkGscore.model$mcmc[[2]][,"intercept"])))
kkgscore.sample <- kintkkgscore[sample.int(nrow(kintkkgscore), 500, replace=FALSE, prob=NULL),]

# roc curve sampling
set.seed(2)
gscore.status.sample <- status.kkgscore[sample.int(nrow(status.kkgscore), 500, replace=FALSE, prob=NULL),]

# make status object for infection prob per time
status.time.KKgscore <- time.steps(status.kkgscore)

t1gs <- as.data.frame(gscore.status.sample [,1:210])
t2gs <- as.data.frame(gscore.status.sample[,211:420])
t3gs <- as.data.frame(gscore.status.sample[,421:630])
t4gs <- as.data.frame(gscore.status.sample[,631:840])

colnames(t1gs) <- gscoreCID
colnames(t2gs) <- gscoreCID
colnames(t3gs) <- gscoreCID
colnames(t4gs) <- gscoreCID

t1gs <- t(t1gs)
t2gs <- t(t2gs)
t3gs <- t(t3gs)
t4gs <- t(t4gs)

t1gs <- as.data.frame(cbind(rownames(t1gs), data.frame(t1gs, row.names=NULL)))
t1gs <- t1gs %>% rename(CID=`rownames(t1gs)`)

t2gs <- as.data.frame(cbind(rownames(t2gs), data.frame(t2gs, row.names=NULL)))
t2gs <- t2gs %>% rename(CID=`rownames(t2gs)`)

t3gs <- as.data.frame(cbind(rownames(t3gs), data.frame(t3gs, row.names=NULL)))
t3gs <- t3gs %>% rename(CID=`rownames(t3gs)`)

t4gs <- as.data.frame(cbind(rownames(t4gs), data.frame(t4gs, row.names=NULL)))
t4gs <- t4gs %>% rename(CID=`rownames(t4gs)`)

#### Gscore ROCR ####

# pre treatment 

predatagscore <- predata %>% dplyr::select(CID, Gscore) %>%
  right_join(t1, by="CID")
predatagscore <- predatagscore[-which(is.na(predatagscore$Gscore)==T),]
gspredictionspreT <- rep(list(predatagscore$Gscore), 500)
gslabelslistpreT <- list()
for(i in 3:ncol(predatagscore)){
  gslabelslistpreT[[i-2]] <- predatagscore[,i]
}

# three weeks 

threeweeksgscore <- threewkdata %>% dplyr::select(CID, Gscore) %>%
  right_join(t2, by="CID")
threeweeksgscore <- threeweeksgscore[-which(is.na(threeweeksgscore$Gscore)==T),]
gspredictionsthreewk <- rep(list(threeweeksgscore$Gscore), 500)
gsabelslistthreewk <- list()
for(i in 3:ncol(threeweeksgscore)){
  gsabelslistthreewk[[i-2]] <- threeweeksgscore[,i]
}

# nine weeks 

nineweeksgscore <- ninewkdata %>% dplyr::select(CID, Gscore) %>%
  right_join(t3, by="CID")
nineweeksgscore <- nineweeksgscore[-which(is.na(nineweeksgscore$Gscore)==T),]
gspredictionsninewk <- rep(list(nineweeksgscore$Gscore), 500)
gslabelslistninewk <- list()
for(i in 3:ncol(nineweeksgscore)){
  gslabelslistninewk[[i-2]] <- nineweeksgscore[,i]
}


# six months 

sixmonthsgscore <- sixmonthdata %>% dplyr::select(CID, Gscore) %>%
  right_join(t4, by="CID")
sixmonthsgscore <- sixmonthsgscore[-which(is.na(sixmonthsgscore$Gscore)==T),]
gspredictionssixmth <- rep(list(sixmonthsgscore$Gscore), 500)
gslabelslistsixmonth <- list()
for(i in 3:ncol(sixmonthsgscore)){
  gslabelslistsixmonth[[i-2]] <- sixmonthsgscore[,i]
}

# plot 
sort(table(cutoffgs3),decreasing=TRUE)[1:6]

manypredgst1 <- prediction(gspredictionspreT, gslabelslistpreT)
auc.perfgst1 <-  performance(manypredgst1, measure = "auc")
aucgs1 <- t(rbindlist(list(auc.perfgst1@y.values)))
quantile(aucgs1,  probs = c(0.025,0.5, 0.975))

sensspec.perfgst1 <- performance(manypredgst1, "sens", "spec")
sens <- list()
spec <- list()
cut <- list()
cutoffgs1 <- vector()
for(i in 1:500){
  sens[[i]] <- slot(sensspec.perfgst1, "y.values")[[i]]
  spec[[i]] <- slot(sensspec.perfgst1, "x.values")[[i]]
  cut[[i]] <- slot(sensspec.perfgst1, "alpha.values")[[i]]
  
  cutoffgs1[i] <- cut[[i]][which.max(sens[[i]] + spec[[i]])]
}

many.roc.gscoret1 <-  performance(manypredgst1, measure = "tpr", x.measure = "fpr")
pdf("pretTrocgscore.pdf") 
plot(many.roc.gscoret1, col="grey82", lty=3)
plot(many.roc.gscoret1, lwd=3, avg="vertical",spread.estimate="boxplot", add=TRUE)
abline(a=0, b= 1)
dev.off()

# three weeks 

manypredgst2 <- prediction(gspredictionsthreewk, gsabelslistthreewk)
auc.perfgst2 <-  performance(manypredgst2, measure = "auc")
aucgs2 <- t(rbindlist(list(auc.perfgst2@y.values)))
quantile(aucgs2,  probs = c(0.025,0.5, 0.975))

sensspec.perfgst2 <- performance(manypredgst2, "sens", "spec")
sens <- list()
spec <- list()
cut <- list()
cutoffgs2 <- vector()
for(i in 1:500){
  sens[[i]] <- slot(sensspec.perfgst2, "y.values")[[i]]
  spec[[i]] <- slot(sensspec.perfgst2, "x.values")[[i]]
  cut[[i]] <- slot(sensspec.perfgst2, "alpha.values")[[i]]
  
  cutoffgs2[i] <- cut[[i]][which.max(sens[[i]] + spec[[i]])]
}

many.roc.gscoret2 <-  performance(manypredgst2, measure = "tpr", x.measure = "fpr")
pdf("threewkrocgscore.pdf") 
plot(many.roc.gscoret2, col="grey82", lty=3)
plot(many.roc.gscoret2, lwd=3, avg="vertical", spread.estimate="boxplot", add=TRUE)
abline(a=0, b= 1)
dev.off()

# nine weeks 
manypredgst3 <- prediction(gspredictionsninewk, gslabelslistninewk)
auc.perfgst3 <-  performance(manypredgst3, measure = "auc")
aucgs3 <- t(rbindlist(list(auc.perfgst3@y.values)))
quantile(aucgs3,  probs = c(0.025,0.5, 0.975))

sensspec.perfgst3 <- performance(manypredgst3, "sens", "spec")
sens <- list()
spec <- list()
cut <- list()
cutoffgs3 <- vector()
for(i in 1:500){
  sens[[i]] <- slot(sensspec.perfgst3, "y.values")[[i]]
  spec[[i]] <- slot(sensspec.perfgst3, "x.values")[[i]]
  cut[[i]] <- slot(sensspec.perfgst3, "alpha.values")[[i]]
  
  cutoffgs3[i] <- cut[[i]][which.max(sens[[i]] + spec[[i]])]
}

many.roc.gscoret3 <-  performance(manypredgst3, measure = "tpr", x.measure = "fpr")
pdf("ninewkrocgscore.pdf") 
plot(many.roc.gscoret3, col="grey82", lty=3)
plot(many.roc.gscoret3, lwd=3, avg="vertical",spread.estimate="boxplot", add=TRUE)
abline(a=0, b= 1)
dev.off()

# six months 
manypredgst4 <- prediction(gspredictionssixmth, gslabelslistsixmonth)
auc.perfgst4 <-  performance(manypredgst4, measure = "auc")
aucgs4 <- t(rbindlist(list(auc.perfgst4@y.values)))
quantile(aucgs4,  probs = c(0.025,0.5, 0.975))

sensspec.perfgst4 <- performance(manypredgst4, "sens", "spec")
sens <- list()
spec <- list()
cut <- list()
cutoffgs4 <- vector()
for(i in 1:500){
  sens[[i]] <- slot(sensspec.perfgst4, "y.values")[[i]]
  spec[[i]] <- slot(sensspec.perfgst4, "x.values")[[i]]
  cut[[i]] <- slot(sensspec.perfgst4, "alpha.values")[[i]]
  
  cutoffgs4[i] <- cut[[i]][which.max(sens[[i]] + spec[[i]])]
}

many.roc.gscoret4 <-  performance(manypredgst4, measure = "tpr", x.measure = "fpr")
pdf("sixmthrocgscore.pdf") 
plot(many.roc.gscoret4, col="grey82", lty=3)
plot(many.roc.gscoret4, lwd=3, avg="vertical", spread.estimate="boxplot", add=TRUE)
abline(a=0, b= 1)
dev.off()

# probability of being infected plot 

pdf("probinfKKgscore.pdf",width=3,height=8)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfrow=c(1,1))

prob.inf.kkgscore <- qplot(value, data = status.time.KKgscore, geom = "histogram", fill=time, binwidth=0.02) +
  scale_fill_manual(values = c("#e7298a", "#1b9e77", "#d95f02", "#7570b3"))+
  geom_text(data=anno, aes(x = xcoord, y = ycoord, label = label))+
  facet_grid(time~.)+
  theme_bw()+coord_cartesian(ylim=c(0, 200), xlim=c(0,1))+theme(legend.position = "none")+
  ylab("Number of Hosts") + xlab("Probability of Being Infected")

prob.inf.kkgscore
dev.off()

# posteriors 

#pdf("k.postsboth.pdf",width=8,height=4)
dev.new()
par(font=2, cex.axis=0.75, lwd=2, mar=c(3,2.5,1.2,0)+0.1,mgp=c(3,0.4,0), bty = "n")
par(mfrow=c(1,2))

plot(klogKKcca$x,klogKKcca$y/max(klogKKcca$y),xlim=c(0,1),ylim=c(0,1),type = "l", ylab = "",xlab="", col="#d95f02")
polygon(c(rev(klogKKcca$x), klogKKcca$x), c(rev(klogKKcca$y/max(klogKKcca$y)), rep(0,length(klogKKcca$y))),
        col = adjustcolor("#d95f02", alpha=.2), border = NA)

mtext("Scaled density",side=2,cex=1,line=1.2)
mtext("logistic function k ",side=1,cex=1,line=1.2)

legend("topright",  "KK & POC-CCA+",
       col=c("#d95f02", "#7570b3"), 
       lty=c(1), bty='n',cex=0.75)
#dev.off()


# cca and gscore intercept plots 

dev.new(height = 6, height = 6)
par(mfrow=c(1,1))
pdf("kposts.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

#plot(klogKKgscore$x,klogKKgscore$y/max(klogKKgscore$y),xlim=c(0,100),ylim=c(0,1),type = "l", ylab = "",xlab="", col="blue")
plot(klogKKcca$x,klogKKcca$y/max(klogKKcca$y),xlim=c(0,100),ylim=c(0,1),type = "l", ylab = "",xlab="", col="red")

mtext("Scaled density",side=2,cex=1,line=1.2)
mtext("logistic function intercept ",side=1,cex=1,line=1.2)
dev.off()
graphics.off()

# logistic curve figure 

# cca
logcurvelistkkcca <- list()
logcurvelistkkgscore <- list()
x=seq(from=1, to=100, by=1)

for(i in 1:100){
  logcurvelistkkcca[[i]] <- as.data.frame(logcurvecca(x,kintkkcca[i,1],kintkkcca[i,2]))
}

for(i in 1:length(logcurvelistkkcca)){
  df <- logcurvelistkkcca[[i]]
  df$x <- as.factor(1:nrow(df))
  logcurvelistkkcca[[i]] <- df
}

for(i in 1:100){
  colnames(logcurvelistkkcca[[i]]) <- c("value", "x")
}

kkccalogcurvesample <- rbindlist(logcurvelistkkcca)

## gscore 
for(i in 1:100){
  logcurvelistkkgscore[[i]] <- as.data.frame(logcurvegscore(x,kintkkgscore[i,1],kintkkgscore[i,2]))
}

for(i in 1:length(logcurvelistkkgscore)){
  df <- logcurvelistkkgscore[[i]]
  df$x <- as.factor(1:nrow(df))
  logcurvelistkkgscore[[i]] <- df
}

for(i in 1:100){
  colnames(logcurvelistkkgscore[[i]]) <- c("value", "x")
}

kkgscorelogcurvesample <- rbindlist(logcurvelistkkgscore)

p <-  c(.025,.5,.975)
p_names <- map_chr(p, ~paste0(.x*100, "%"))

p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

kkccaquantiles <- kkccalogcurvesample %>% group_by(x) %>%
  summarise_at(vars(value),funs(!!!p_funs))

kkgscorequantiles <- kkgscorelogcurvesample %>% group_by(x) %>%
  summarise_at(vars(value),funs(!!!p_funs))

kkccalog <- logcurve(x,mean(klogKKcca$x), mean(interceptKKcca$x))
kkgscorelog <- logcurve(x,mean(klogKKgscore$x), mean(interceptKKgscore$x))

x=seq(from=0, to=100, by =1)
polygonx <- c(kkccaquantiles$x, rev(kkccaquantiles$x))
polygony <- c(kkccaquantiles$`2.5%`, rev(kkccaquantiles$`97.5%`))

polygonxgscore <- c(kkgscorequantiles$x, rev(kkgscorequantiles$x))
polygonygscore <- c(kkgscorequantiles$`2.5%`, rev(kkgscorequantiles$`97.5%`))


par(mfrow=c(1,2))
plot(NA,NA, xlim = c(0,100), ylim = c(0,9), xlab="Kato-Katz Count", ylab="f(x)", cex.axis=1)
lines(kkgscorequantiles$`2.5%`, col="#8da0cb", lwd=1) 
polygon(x=polygonxgscore, y=polygonygscore, col=adjustcolor("grey", alpha.f=0.4), border=NA)
axis(2, seq(0,9,1))
plot(NA,NA, xlim = c(0,100), ylim = c(0,4), xlab="Kato-Katz Count", ylab="f(x)", cex.axis=1)
lines(kkccaquantiles$`50%`, col="#fc8d62", lwd=2, lty = 1)
polygon(x=polygonx, y=polygony, col=adjustcolor("grey", alpha.f=0.4), border=NA)


dev.copy(pdf, "logcurve.confint.28.8.20.pdf", width = 12, height = 6)
dev.off()
graphics.off()

#### log curve with data and noise ####
logvecdatacca <- vector()
logvecdatags <- vector()

mergekkcca$logvecdatacca <- logcurvecca(mergekkcca$mean.epg,mean(kkcca.sample[,1]),mean(kkcca.sample[,2]))
mergekkcca$kkgscorelog <- logcurvegscore(mergekkcca$mean.epg,mean(kkgscore.sample[,1]), mean(kkgscore.sample[,1]))

mergekkcca$cca[which(mergekkcca$Poppy=="+++")]<-4
mergekkcca$cca[which(mergekkcca$Poppy== "++")]<-3
mergekkcca$cca[which(mergekkcca$Poppy=="+")]<- 2
mergekkcca$cca[which(mergekkcca$Poppy=="Trace")]<- 1
mergekkcca$cca[which(mergekkcca$Poppy=="Negative (0)")]<- 0
mergekkcca$sdmin <- mergekkcca$logvecdatacca-1.758821
mergekkcca$sdmax <- mergekkcca$logvecdatacca+1.758821
mergekkcca$sdmax[which(mergekkcca$sdmax>4)] <- 4

mergekkcca$gssdmin <- mergekkcca$kkgscorelog-1.045756
mergekkcca$gssdmax <- (mergekkcca$kkgscorelog+1.045756)+1
mergekkcca$gssdmax[which(mergekkcca$gssdmax>10)] <- 10

ccalog <- ggplot(data=mergekkcca)+
  geom_point(aes(x=mean.epg, y=cca), colour="darkred", shape=1)+
  geom_ribbon(aes(x=mean.epg, ymin=sdmin, ymax=sdmax), alpha=0.2, fill="#fc8d62")+
  geom_line(aes(x=mean.epg ,y=logvecdatacca), colour="darkred", lwd=1)+
  geom_vline(xintercept = 99, lty=2, colour="lightgrey")+
  geom_vline(xintercept = 399, lty=2, colour="lightgrey")+
  xlim(0,450)+scale_y_continuous(breaks=c(0:4))+
  theme_classic()+
  ylab("CCA+")+xlab("Mean Eggs Per Gram")

ccalogzoom <- ggplot(data=mergekkcca)+
  geom_point(aes(x=mean.epg, y=cca), colour="darkred", shape=1)+
  geom_ribbon(aes(x=mean.epg, ymin=sdmin, ymax=sdmax), alpha=0.2, fill="#fc8d62")+
  geom_line(aes(x=mean.epg ,y=logvecdatacca), colour="darkred", lwd=1)+
  #geom_vline(xintercept = 99, lty=2, colour="lightgrey")+
  #geom_vline(xintercept = 399, lty=2, colour="lightgrey")+
  xlim(0,20)+scale_y_continuous(breaks=c(0:4))+
  theme_classic()+
  ylab("CCA+")+xlab("Mean Eggs Per Gram")

gslog <- ggplot(data=mergekkcca)+
  geom_point(aes(x=mean.epg, y=gscore), colour="#8da0cb", shape=1)+
  geom_ribbon(aes(x=mean.epg, ymin=gssdmin+1, ymax=gssdmax), alpha=0.3, fill="#8da0cb")+
  geom_line(aes(x=mean.epg ,y=kkgscorelog+1), colour="#8da0cb", lwd=1)+
  geom_vline(xintercept = 99, lty=2, colour="lightgrey")+
  geom_vline(xintercept = 399, lty=2, colour="lightgrey")+
  xlim(0,450)+scale_y_continuous(breaks=c(1:10))+
  theme_classic()+
  ylab("G-Score")+xlab("Mean Eggs Per Gram")

gslogzoom <- ggplot(data=mergekkcca)+
  geom_point(aes(x=mean.epg, y=gscore), colour="#8da0cb", shape=1)+
  geom_ribbon(aes(x=mean.epg, ymin=gssdmin+1, ymax=gssdmax), alpha=0.3, fill="#8da0cb")+
  geom_line(aes(x=mean.epg ,y=kkgscorelog+1), colour="#8da0cb", lwd=1)+
  #geom_vline(xintercept = 99, lty=2, colour="lightgrey")+
  #geom_vline(xintercept = 399, lty=2, colour="lightgrey")+
  xlim(0,20)+scale_y_continuous(breaks=c(1:10))+
  theme_classic()+
  ylab("G-Score")+xlab("Mean Eggs Per Gram")

logpanel <- ggarrange(ccalog, gslog,ccalogzoom,gslogzoom, labels=c("A","B", "C", "D"), nrow=2, ncol=2)+ggsave("logpanel.pdf")

logcurve <- logcurvegscore(x=seq(from=1, to=400, by=1), k=mean(klogKKgscore),intercept=mean(interceptKKgscore))

#### just KK model ####

katokatz <- CalModKK(dim(kk)[1], dim(kk)[2], dim(kk)[3], kk)

prevKK <- density(c(katokatz$mcmc[[1]][,"prev"], katokatz$mcmc[[2]][,"prev"]))
model.output <- as.data.frame(as.matrix(as.mcmc(katokatz)))
status <- model.output[,3:842]
status.time.KK <- time.steps(status)
# roc curve sampling
set.seed(2)
kk.status.sample <- status[sample.int(nrow(status), 500, replace=FALSE, prob=NULL),]

t1kk <- as.data.frame(kk.status.sample [,1:210])
t2kk <- as.data.frame(kk.status.sample[,211:420])
t3kk <- as.data.frame(kk.status.sample[,421:630])
t4kk <- as.data.frame(kk.status.sample[,631:840])

colnames(t1kk) <- CID
colnames(t2kk) <- CID
colnames(t3kk) <- CID
colnames(t4kk) <- CID

t1kk <- t(t1kk)
t2kk <- t(t2kk)
t3kk <- t(t3kk)
t4kk <- t(t4kk)

t1kk <- as.data.frame(cbind(rownames(t1kk), data.frame(t1kk, row.names=NULL)))
t1kk <- t1kk %>% rename(CID=`rownames(t1kk)`)

t2kk <- as.data.frame(cbind(rownames(t2kk), data.frame(t2kk, row.names=NULL)))
t2kk <- t2kk %>% rename(CID=`rownames(t2kk)`)

t3kk <- as.data.frame(cbind(rownames(t3kk), data.frame(t3kk, row.names=NULL)))
t3kk <- t3kk %>% rename(CID=`rownames(t3kk)`)

t4kk <- as.data.frame(cbind(rownames(t4kk), data.frame(t4kk, row.names=NULL)))
t4kk <- t4kk %>% rename(CID=`rownames(t4kk)`)


# This is a df for the geom_text with the co-ordinates of the labels 
anno <- data.frame(c("Baseline", "ThreeWeeks", "NineWeeks", "SixMonths"))
colnames(anno) <- "time"
anno$xcoord <- 0
anno$ycoord <-200
anno$label <- as.factor(c("A", "B", "C", "D"))

pdf("probinfKK.pdf",width=3,height=8)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfrow=c(1,1))

prob.inf.kk <- qplot(value, data = status.time.KK, geom = "histogram", fill=time, binwidth=0.02) +
  scale_fill_manual(values = c("#e7298a", "#1b9e77", "#d95f02", "#7570b3"))+
  geom_text(data=anno, aes(x = xcoord, y = ycoord, label = label))+
  facet_grid(time~.)+
  theme_bw()+coord_cartesian(ylim=c(0, 200), xlim=c(0,1))+theme(legend.position = "none")+
  ylab("Number of Hosts") + xlab("Probability of Being Infected")

prob.inf.kk
dev.off()

# all the models prob inf panel 
status.time.KKcca <- status.time.KKcca %>% mutate(diagnostic="Plus(+)")
status.time.KKgscore <- status.time.KKgscore %>% mutate(diagnostic="G-Score")
status.time.KK <- status.time.KK %>% mutate(diagnostic="Kato-Katz")

probinfdf <- bind_rows(status.time.KKcca,status.time.KKgscore, status.time.KK )%>%
  mutate_if(is.character, as.factor)%>%
  mutate(diagnostic=factor(diagnostic, levels=c("Kato-Katz", "Plus(+)", "G-Score")))

cols=c("#66c2a5", "#fc8d62", "#8da0cb")
prob.inf.all.plot <- ggplot()+
  geom_histogram(data=probinfdf, aes(x=value, fill=diagnostic), bins=50)+
  facet_grid(time~diagnostic)+
  scale_fill_manual(values=cols)+
  ylab("Number of Children") + xlab("Probability of Being Infected")+
  theme_bw()+
  ggsave("prob.inf.all.pdf")

#### raw data plots ####

# restructure for categorical
raw.ccadata <- ccadata
raw.ccadata$Poppy <- as.character(raw.ccadata$Poppy)
raw.ccadata$Poppy[which(raw.ccadata$Poppy=="T")] <- "Trace"
raw.ccadata$Poppy[which(raw.ccadata$Poppy==3)]<-"+++"
raw.ccadata$Poppy[which(raw.ccadata$Poppy==2)]<- "++"
raw.ccadata$Poppy[which(raw.ccadata$Poppy==1)]<- "+"
raw.ccadata$Poppy[which(raw.ccadata$Poppy==0.5)]<- "Trace"
raw.ccadata$Poppy[which(raw.ccadata$Poppy==0)]<- "Negative (0)"
raw.ccadata$Poppy <- factor(raw.ccadata$Poppy)
raw.ccadata <- raw.ccadata %>% filter(Poppy!="-")
raw.ccadata$Poppy <- factor(raw.ccadata$Poppy)
raw.ccadata <- raw.ccadata %>% filter(Poppy!="")
raw.ccadata$Poppy <- factor(raw.ccadata$Poppy, levels=c("Negative (0)", "Trace", "+", "++","+++"))
raw.ccadata$gscore <- as.integer(raw.ccadata$gscore)

raw.ccadata$dateN <- NA
raw.ccadata$dateN[which(raw.ccadata$date_on_tube=="25/09/2017" | raw.ccadata$date_on_tube=="26/09/2017" | raw.ccadata$date_on_tube=="27/09/2017"
                         | raw.ccadata$date_on_tube=="28/09/2017" | raw.ccadata$date_on_tube=="29/09/2017" | raw.ccadata$date_on_tube=="02/10/2017")] <- "Pre-T"

raw.ccadata$dateN[which(raw.ccadata$date_on_tube=="23/10/2017" | raw.ccadata$date_on_tube=="24/10/2017" | raw.ccadata$date_on_tube=="25/10/2017"
                         | raw.ccadata$date_on_tube=="26/10/2017" | raw.ccadata$date_on_tube=="27/10/2017")] <- "3 weeks"

raw.ccadata$dateN[which(raw.ccadata$date_on_tube=="04/12/2017" | raw.ccadata$date_on_tube=="05/12/2017")] <- "9 weeks"

raw.ccadata$dateN[which(raw.ccadata$date_on_tube=="01/03/2018" | raw.ccadata$date_on_tube=="05/03/2018" | raw.ccadata$date_on_tube=="06/03/2018"
                         | raw.ccadata$date_on_tube=="07/03/2018" | raw.ccadata$date_on_tube=="08/03/2018" | raw.ccadata$date_on_tube=="09/03/2018")] <- "6 months"

ccadates <- raw.ccadata[,c("CID", "Poppy", "gscore", "dateN")]
ccadates$dateN <- as.factor(ccadates$dateN)

# sort out KK data for plot 

kkxgscore.data <- kkdata
kkxgscore.data$dateN <- NA
kkxgscore.data$dateN[which(kkxgscore.data$date=="25/09/2017" | kkxgscore.data$date=="26/09/2017" | kkxgscore.data$date=="27/09/2017"
                           | kkxgscore.data$date=="28/09/2017" | kkxgscore.data$date=="29/09/2017" | kkxgscore.data$date=="02/10/2017")] <- "Pre-T"

kkxgscore.data$dateN[which(kkxgscore.data$date=="23/10/2017" | kkxgscore.data$date=="24/10/2017" | kkxgscore.data$date=="25/10/2017"
                           | kkxgscore.data$date=="26/10/2017" | kkxgscore.data$date=="27/10/2017")] <- "3 weeks"

kkxgscore.data$dateN[which(kkxgscore.data$date=="04/12/2017" | kkxgscore.data$date=="05/12/2017")] <- "9 weeks"

kkxgscore.data$dateN[which(kkxgscore.data$date=="01/03/2018" | kkxgscore.data$date=="05/03/2018" | kkxgscore.data$date=="06/03/2018"
                           | kkxgscore.data$date=="07/03/2018" | kkxgscore.data$date=="08/03/2018" | kkxgscore.data$date=="09/03/2018")] <- "6 months"

kkxgscore.data <- melt(kkxgscore.data,id.vars = c("child_id", "date", "dateN"), measure.vars = c("Sm_A", "Sm_B"))
kkxgscore.data$dateN <- as.factor(kkxgscore.data$dateN)
kkxgscore.data <- kkxgscore.data[-which(is.na(kkxgscore.data$value)==T),]

kkxgscore.data$child_id <- as.factor(kkxgscore.data$child_id)

kkxgscore.data <- kkxgscore.data %>% 
  group_by(dateN, child_id)%>%
  summarise(mean_count=mean(value, na.rm=T))

colnames(kkxgscore.data) <- c("dateN", "CID", "mean_count")

# merge data 

mergekkcca <- merge(ccadates, kkxgscore.data, by=c("CID", "dateN"))
mergekkcca$gscore <- as.factor(mergekkcca$gscore)
mergekkcca$dateN <- factor(mergekkcca$dateN, levels=c("Pre-T", "3 weeks", "9 weeks", "6 months"))
mergekkcca$logcount <- log(mergekkcca$mean_count*24)
mergekkcca[mergekkcca=="-Inf"] <- 0
mergekkcca <- mergekkcca %>% mutate(mean.epg=mean_count*24, intcat=NA)
for(i in 1:nrow(mergekkcca)){
  if(mergekkcca[i,ncol(mergekkcca)-1]==0){
    mergekkcca[i,ncol(mergekkcca)] <- "zero_count"
  } else if(mergekkcca[i,ncol(mergekkcca)-1]>0&&mergekkcca[i,ncol(mergekkcca)-1]<=99){
    mergekkcca[i,ncol(mergekkcca)] <- "low"
  } else if(mergekkcca[i,ncol(mergekkcca)-1]>=100&&mergekkcca[i,ncol(mergekkcca)-1]<=399){
    mergekkcca[i,ncol(mergekkcca)] <- "moderate"
  } else {
    mergekkcca[i,ncol(mergekkcca)] <- "high"
  }
}
mergekkcca$intcat <- as.factor(mergekkcca$intcat)  
mergekkcca$intcat <- factor(mergekkcca$intcat, levels = c("zero_count", "low", "moderate", "high"))

#mergekkcca$Poppy <- as.numeric(mergekkcca$Poppy)
mergekkcca$gscore <- as.numeric(mergekkcca$gscore)

corr.diag.plot <- ggplot(data=mergekkcca, aes(x=gscore, y=Poppy))+
  geom_jitter(aes(colour=intcat), size=2)+
  facet_grid(dateN~.)+
  labs(colour="Intensity Categories", y="Plus (+) Score", x="G-Score")+
  scale_colour_viridis(discrete=TRUE) +
  theme_bw()+
  scale_x_continuous(breaks=c(1:10,1))+
  ggsave("corr.diag.plot.raw.pdf")

#### prob infected by raw kk count by CCA value ####

kkpre$CID <- as.factor(kkpre$CID)

probinfccalist <- split(status.time.KKcca, status.time.KKcca$time)
for(i in 1:length(probinfccalist)){
  colnames(probinfccalist[[i]]) <- c("time", "CID", "prob.inf")
}

#### CCA +
probinfccalist[[1]] <- kkpre %>% dplyr::select( CID, meancount)%>%
  right_join(probinfccalist[[1]], by="CID")

probinfccalist[[1]] <- ccapre %>% dplyr::select(CID, CCA)%>%
  right_join(probinfccalist[[1]], by="CID")

probinfccalist[[2]] <- kkthreewk %>% dplyr::select( CID, meancount)%>%
  right_join(probinfccalist[[2]], by="CID")

probinfccalist[[2]] <- ccathreewk %>% dplyr::select(CID, CCA)%>%
  right_join(probinfccalist[[2]], by="CID")

probinfccalist[[3]] <- kkninewk %>% dplyr::select( CID, meancount)%>%
  right_join(probinfccalist[[3]], by="CID")

probinfccalist[[3]] <- ccaninewk %>% dplyr::select(CID, CCA)%>%
  right_join(probinfccalist[[3]], by="CID")

probinfccalist[[4]] <- kksixmth %>% dplyr::select( CID, meancount)%>%
  right_join(probinfccalist[[4]], by="CID")

probinfccalist[[4]] <- ccasixmth %>% dplyr::select(CID, CCA)%>%
  right_join(probinfccalist[[4]], by="CID")

probinfcca <- rbindlist(probinfccalist)%>%
  mutate_if(is.character, as.factor)
probinfcca$CCA <- as.factor(probinfcca$CCA)
probinfcca <- probinfcca[-which(is.na(probinfcca$CCA)==T),]

col=c("#edf8fb","#b3cde3","#8c96c6","#8856a7", "#810f7c")
col=c("#fdd0a2","#fd8d3c","#f16913","#d94801","#8c2d04")
ggplot(data=probinfcca)+
  geom_jitter(aes(x=prob.inf, y=meancount, colour=CCA), size=3)+
  facet_grid(time~.)+
  labs(colour="CCA+ Score", y="Mean Kato-Katz Count", x="Probability Infected")+
  scale_colour_manual(values=col)+
  theme_bw()


probinfcca$`eggs per gram` <- probinfcca$meancount*24 
probinfcca <- probinfcca %>% mutate(infection_category = case_when(`eggs per gram`==0 | `eggs per gram`==0.0 ~ "zero count", 
                                        `eggs per gram`>0 & `eggs per gram` <=99~"light infection", 
                                        `eggs per gram`>=100&`eggs per gram`<=399 ~"moderate infection", 
                                        TRUE~"heavy infection"))%>%
                                      mutate_if(is.character, as.factor)%>%
                                      mutate(infection_category = factor(infection_category, levels=c("zero count","light infection", "moderate infection", "heavy infection" )))

col=c("#fdd0a2","#fd8d3c","#f16913","#d94801","#8c2d04")
ggplot(data=probinfcca )+
  geom_count(aes(x=infection_category, y=CCA, color = prob.inf, size = ..n..)) +
  guides(color = 'legend')+
  #geom_jitter(aes(x=infection_category, y=CCA, colour=prob.inf), size=2)+
  facet_grid(time~.)+
  labs(colour="Probability Infected", y="CCA +", x="WHO Infection Intensity Category")+
  theme_bw()


#### CCA G-Score

probinfgscorelist <- split(status.time.KKgscore, status.time.KKgscore$time)

probinfgscorelist <- split(status.time.KKcca, status.time.KKcca$time)
for(i in 1:length(probinfgscorelist)){
  colnames(probinfgscorelist[[i]]) <- c("time", "CID", "prob.inf")
}

probinfgscorelist[[1]] <- kkpre %>% dplyr::select( CID, meancount)%>%
  right_join(probinfgscorelist[[1]], by="CID")

probinfgscorelist[[1]] <- gscorepre %>% dplyr::select(CID, Gscore)%>%
  right_join(probinfgscorelist[[1]], by="CID")

probinfgscorelist[[2]] <- kkthreewk %>% dplyr::select( CID, meancount)%>%
  right_join(probinfgscorelist[[2]], by="CID")

probinfgscorelist[[2]] <- gscorethreewk %>% dplyr::select(CID, Gscore)%>%
  right_join(probinfgscorelist[[2]], by="CID")

probinfgscorelist[[3]] <- kkninewk %>% dplyr::select( CID, meancount)%>%
  right_join(probinfgscorelist[[3]], by="CID")

probinfgscorelist[[3]] <- gscoreninewk %>% dplyr::select(CID, Gscore)%>%
  right_join(probinfgscorelist[[3]], by="CID")

probinfgscorelist[[4]] <- kksixmth %>% dplyr::select( CID, meancount)%>%
  right_join(probinfgscorelist[[4]], by="CID")

probinfgscorelist[[4]] <- gscoresixmth %>% dplyr::select(CID, Gscore)%>%
  right_join(probinfgscorelist[[4]], by="CID")

probinfgscore <- rbindlist(probinfgscorelist)


#### prob inf x cca score x WHO cats 
probinfmodlist <- split(probinfdf, probinfdf$diagnostic)
probinfmodlist <- lapply(probinfmodlist, function(x) split(x, x$time))

for(i in 1:3){
    for(j in 1:4){
      probinfmodlist[[i]][[j]] <- data.frame(probinfmodlist[[i]][[j]])
      colnames(probinfmodlist[[i]][[j]]) <- c("time", "CID", "prob.inf", "diagnostic")
    }
}

probinfmodlist[[1]][[1]] <- kkpre %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[1]][[1]], by="CID")

probinfmodlist[[1]][[2]] <- kkthreewk %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[1]][[2]], by="CID")

probinfmodlist[[1]][[3]] <- kkninewk %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[1]][[3]], by="CID")

probinfmodlist[[1]][[4]] <- kksixmth %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[1]][[4]], by="CID")

probinfKK2 <- bind_rows(probinfmodlist[[1]][[1]], probinfmodlist[[1]][[2]], probinfmodlist[[1]][[3]], probinfmodlist[[1]][[4]])



probinfmodlist[[2]][[1]] <- ccapre %>% dplyr::select(CID, CCA) %>%
  right_join(probinfmodlist[[2]][[1]], by="CID")

probinfmodlist[[2]][[1]] <- kkpre %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[2]][[1]], by="CID")

probinfmodlist[[2]][[2]] <- ccathreewk %>% dplyr::select(CID, CCA) %>%
  right_join(probinfmodlist[[2]][[2]], by="CID")

probinfmodlist[[2]][[2]] <- kkthreewk %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[2]][[2]], by="CID")

probinfmodlist[[2]][[3]] <- ccaninewk %>% dplyr::select(CID, CCA) %>%
  right_join(probinfmodlist[[2]][[3]], by="CID")

probinfmodlist[[2]][[3]] <- kkninewk %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[2]][[3]], by="CID")

probinfmodlist[[2]][[4]] <- ccasixmth %>% dplyr::select(CID, CCA) %>%
  right_join(probinfmodlist[[2]][[4]], by="CID")

probinfmodlist[[2]][[4]] <- kksixmth %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[2]][[4]], by="CID")

probinfcca2 <- bind_rows(probinfmodlist[[2]][[1]], probinfmodlist[[2]][[2]], probinfmodlist[[2]][[3]], probinfmodlist[[2]][[4]])
probinfcca2$`eggs per gram` <- probinfcca2$meancount*24

probinfcca2 <- probinfcca2 %>% mutate(infection_category = case_when(`eggs per gram`==0 | `eggs per gram`==0.0 ~ "zero count", 
                                                                     `eggs per gram`>0 & `eggs per gram` <=99~"light infection", 
                                                                     `eggs per gram`>=100&`eggs per gram`<=399 ~"moderate infection", 
                                                                     TRUE~"heavy infection"))%>%
  mutate_if(is.character, as.factor)%>%
  mutate(infection_category = factor(infection_category, levels=c("zero count","light infection", "moderate infection", "heavy infection" )),
         CCA = factor(CCA, levels=c("0", "1", "2", "3", "4"), labels=c("Negative", "Trace", "+", "++", "+++")))

probinfcca2$CCA <- as.factor(probinfcca2$CCA)
probinfcca2 <- probinfcca2[-which(is.na(probinfcca2$CCA)==T),]

col=c("#fdd0a2","#fd8d3c","#f16913","#d94801","#8c2d04")
# too light col=c("#9970ab","#c2a5cf","#d9f0d3","#a6dba0","#5aae61")
ggplot(data=probinfcca2)+
  geom_jitter(aes(x=prob.inf, y=CCA, colour=CCA),size=2.4, shape=1, stroke = 1.2)+
  facet_grid(time~infection_category)+
  labs(colour="CCA+ Score", y="CCA+ Score", x="Probablity Infected")+
  scale_colour_manual(values=col)+
  theme_bw()+
  ggsave("cca.calibration.facet.plot2910.pdf")

probinfmodlist[[3]][[1]] <- gscorepre %>% dplyr::select(CID, Gscore) %>%
  right_join(probinfmodlist[[3]][[1]], by="CID")

probinfmodlist[[3]][[1]] <- kkpre %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[3]][[1]], by="CID")

probinfmodlist[[3]][[1]] <- ccapre %>% dplyr::select(CID, CCA) %>%
  right_join(probinfmodlist[[3]][[1]], by="CID")

probinfmodlist[[3]][[2]] <- gscorethreewk %>% dplyr::select(CID, Gscore) %>%
  right_join(probinfmodlist[[3]][[2]], by="CID")

probinfmodlist[[3]][[2]] <- kkthreewk %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[3]][[2]], by="CID")

probinfmodlist[[3]][[2]] <- ccapre %>% dplyr::select(CID, CCA) %>%
  right_join(probinfmodlist[[3]][[2]], by="CID")

probinfmodlist[[3]][[3]] <- gscoreninewk %>% dplyr::select(CID, Gscore) %>%
  right_join(probinfmodlist[[3]][[3]], by="CID")

probinfmodlist[[3]][[3]] <- kkninewk %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[3]][[3]], by="CID")

probinfmodlist[[3]][[3]] <- ccapre %>% dplyr::select(CID, CCA) %>%
  right_join(probinfmodlist[[3]][[3]], by="CID")

probinfmodlist[[3]][[4]] <- gscoresixmth %>% dplyr::select(CID, Gscore) %>%
  right_join(probinfmodlist[[3]][[4]], by="CID")

probinfmodlist[[3]][[4]] <- kksixmth %>% dplyr::select(CID, meancount) %>%
  right_join(probinfmodlist[[3]][[4]], by="CID")

probinfmodlist[[3]][[4]] <- ccapre %>% dplyr::select(CID, CCA) %>%
  right_join(probinfmodlist[[3]][[4]], by="CID")

probinfgscore2 <- bind_rows(probinfmodlist[[3]][[1]], probinfmodlist[[3]][[2]], probinfmodlist[[3]][[3]], probinfmodlist[[3]][[4]])

probinfgscore2$`eggs per gram` <- probinfgscore2$meancount*24

probinfgscore2 <- probinfgscore2 %>% mutate(infection_category = case_when(`eggs per gram`==0 | `eggs per gram`==0.0 ~ "zero count", 
                                                                     `eggs per gram`>0 & `eggs per gram` <=99~"light infection", 
                                                                     `eggs per gram`>=100&`eggs per gram`<=399 ~"moderate infection", 
                                                                     TRUE~"heavy infection"))%>%
  mutate_if(is.character, as.factor)%>%
  mutate(infection_category = factor(infection_category, levels=c("zero count","light infection", "moderate infection", "heavy infection" )),
         Gscore = factor(Gscore, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), labels=c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10")), 
         CCA = factor(CCA, levels=c("0", "1", "2", "3", "4"), labels=c("Negative", "Trace", "+", "++", "+++")))

probinfgscore2$Gscore <- as.factor(probinfgscore2$Gscore)
probinfgscore2$CCA <- as.factor(probinfgscore2$CCA) 
probinfgscore2 <- probinfgscore2[-which(is.na(probinfgscore2$Gscore)==T),]

probinfgscore3 <- probinfgscore2[-which(is.na(probinfgscore2$CCA)==T),]

col=c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#c7eae5","#80cdc1","#35978f","#01665e","#003c30")
ggplot(data=probinfgscore2)+
  geom_jitter(aes(x=prob.inf, y=Gscore, colour=Gscore), size=2.4, shape=1, stroke = 1.2)+
  facet_grid(time~infection_category)+
  labs(colour="G-Score", y="G-Score", x="Probability Infected")+
  #scale_color_viridis(option="cividis",discrete=TRUE) +
  scale_colour_manual(values=col)+
  theme_bw()+
  ggsave("gscore.calibration.facet.plot2910.pdf")

#### combined plot with cca x gscore and prob infected with facet time and inf cat

ggplot(data=probinfgscore3)+
  geom_jitter(aes(x=CCA, y=Gscore, colour=prob.inf), shape = 1, stroke=1)+
  facet_grid(time~infection_category)+
  labs(colour="Probability Infeted", y="G-Score", x="CCA+")+
  #scale_color_viridis(option="cividis",discrete=TRUE) +
  #scale_colour_manual(values=col)+
  theme_bw()+scale_colour_gradientn(colours=rainbow(5))+
  ggsave("probinfcorrplot.pdf")


mergeprobs <- mergekkcca
mergeprobs <- mergeprobs %>% rename(time=dateN)%>%
  mutate(time=factor(time, levels=c("Pre-T", "3 weeks", "9 weeks", "6 months"), labels = c("Baseline", "ThreeWeeks", "NineWeeks", "SixMonths")))
status.time.KKcca <- rename(status.time.KKcca, CID=variable)
mergeprobs <- status.time.KKcca %>% dplyr::select(time, CID, value)%>%right_join(mergeprobs, by=c("time","CID" ))%>%
  rename(PI.plusmod=value)

status.time.KKgscore<- rename(status.time.KKgscore, CID=variable)

mergeprobs <- status.time.KKgscore %>% dplyr::select(time, CID, value)%>%right_join(mergeprobs, by=c("time","CID" ))%>%
  rename(PI.gscoremod=value)
mergeprobs$Poppy <- as.factor(mergeprobs$Poppy)
mergeprobs$gscore <- as.factor(mergeprobs$gscore)

mergeprobs <- melt(mergeprobs, measure.vars = c("PI.gscoremod", "PI.plusmod"), id.vars = c("time", "CID", "Poppy", "gscore", "mean_count", "logcount", "mean.epg", "intcat"))

mergeprobsgscore <- mergeprobs %>% filter(variable=="PI.gscoremod")
mergeprobscca <- mergeprobs %>% filter(variable=="PI.plusmod")

# set seed to make jitter reproducible 
#set.seed(1); mergeprobsgscore %>% ggplot()+
gscoreprobsplot <-   ggplot(data=mergeprobsgscore)+
  geom_jitter(aes(x=intcat, y=gscore, colour=value), shape = 1, stroke=1)+
  facet_grid(time~.)+
  labs(colour="Probability Infected", y="G-Score", x="WHO Infection Category")+
  theme_bw()+scale_colour_gradientn(colours=rainbow(5))

#set.seed(1); mergeprobscca %>% ggplot()+
ccaprobsplot <- ggplot(data=mergeprobscca)+
  geom_jitter(aes(x=intcat, y=Poppy, colour=value), shape = 1, stroke=1)+
  facet_grid(time~.)+
  labs(colour="Probability Infected", y="CCA +", x="WHO Infection Category")+
  theme_bw()+scale_colour_gradientn(colours=rainbow(5))

 probs.plots <- ggarrange(gscoreprobsplot, ccaprobsplot, labels="AUTO", common.legend = T)+ggsave("probplots.pdf")

 countscats <- mergeprobs%>% group_by(intcat, Poppy, time)%>%tally()

ggplot(data=countscats,aes(x=intcat, y=n, group=Poppy) )+
  geom_col(aes(fill=Poppy), position ="dodge")+
  geom_text(aes(label=n, y=n+0.05), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_grid(time~.)+
  labs(fill="CCA+ Scores", y="Number of Children +", x="WHO Infection Category")+
  theme_bw()+
  ggsave("dist.ccascoresxinfcat.pdf")

countscatsgs <- mergeprobs%>% group_by(intcat, gscore, time)%>%tally()

ggplot(data=countscatsgs,aes(x=intcat, y=n, group=gscore) )+
  geom_col(aes(fill=gscore), position ="dodge")+
  geom_text(aes(label=n, y=n+0.05), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_grid(time~.)+
  labs(fill="G-Scores", y="Number of Children +", x="WHO Infection Category")+
  theme_bw()+
  ggsave("dist.gsscoresxinfcat.pdf")

# prevalence plot #

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="Prevalence", ylab="Scaled Density", cex.axis=1)
lines(prevKK$x,prevKK$y/max(prevKK$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(prevKK$x), prevKK$x), c(rev(prevKK$y/max(prevKK$y)), rep(0,length(prevKK$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(prevKKGScore$x,prevKKGScore$y/max(prevKKGScore$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(prevKKGScore$x), prevKKGScore$x), c(rev(prevKKGScore$y/max(prevKKGScore$y)), rep(0,length(prevKKGScore$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(prevKKGScore2$x,prevKKGScore2$y/max(prevKKGScore2$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(prevKKGScore2$x), prevKKGScore2$x), c(rev(prevKKGScore2$y/max(prevKKGScore2$y)), rep(0,length(prevKKGScore2$y))),
        col = adjustcolor("black", alpha=.3))
lines(prevKKcca$x,prevKKcca$y/max(prevKKcca$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(prevKKcca$x), prevKKcca$x), c(rev(prevKKcca$y/max(prevKKcca$y)), rep(0,length(prevKKcca$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topleft", c("Kato-Katz", "KK & G-Score", "KK & Plus (+)"),
       col=c("#1b9e77","#d95f02", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "prevalence.20.7.20.pdf", height = 6, width = 6)
dev.off()
graphics.off()


#### cumulative percentage plots 

ccacumlist <- list()
for(i in 2:ncol(t1cca)){
  ccacumlist[[i-1]] <- as.data.frame(t1cca[,i])
  colnames(ccacumlist[[i-1]]) <-"status" 
  noinf <- filter(ccacumlist[[i-1]],status==0)
  noinf$rank <- 0
  ccacumlist[[i-1]] <- ccacumlist[[i-1]][ccacumlist[[i-1]]$status!=0,]
  ccacumlist[[i-1]] <- as.data.frame(ccacumlist[[i-1]])
  ccacumlist[[i-1]]$rank <- 1:nrow(ccacumlist[[i-1]])
  ccacumlist[[i-1]] <- bind_rows(noinf, ccacumlist[[i-1]])
  ccacumlist[[i-1]]$cumtot <- (ccacumlist[[i-1]]$rank/210)*100
}

ccacumarray <- array(as.numeric(unlist(ccacumlist)), dim=c(210,4,500))

ccacumquants <- as.data.frame(matrix(nrow=210, ncol = 3))
#take the quantiles across the sampled 500
for(j in 1:210){
  ccacumquants[j,] <- quantile(ccacumarray[j,4,1:500], probs=c(0.025, 0.5, 0.975))
}

gscorecumlist <- list()
for(i in 2:ncol(t1gs)){
  gscorecumlist[[i-1]] <- as.data.frame(t1gs[,i])
  colnames(gscorecumlist[[i-1]]) <-"status" 
  noinf <- filter(gscorecumlist[[i-1]],status==0)
  noinf$rank <- 0
  gscorecumlist[[i-1]] <- gscorecumlist[[i-1]][gscorecumlist[[i-1]]$status!=0,]
  gscorecumlist[[i-1]] <- as.data.frame(gscorecumlist[[i-1]])
  gscorecumlist[[i-1]]$rank <- 1:nrow(gscorecumlist[[i-1]])
  gscorecumlist[[i-1]] <- bind_rows(noinf, gscorecumlist[[i-1]])
  gscorecumlist[[i-1]]$cumtot <- (gscorecumlist[[i-1]]$rank/210)*100
}

gscumarray <- array(as.numeric(unlist(gscorecumlist)), dim=c(210,4,500))
gscorecumquants <- as.data.frame(matrix(nrow=210, ncol = 3))
#take the quantiles across the sampled 500
for(j in 1:210){
  gscorecumquants[j,] <- quantile(gscumarray[j,4,1:500], probs=c(0.025, 0.5, 0.975))
}

kkcumlist <- list()
for(i in 2:ncol(t1kk)){
  kkcumlist[[i-1]] <- as.data.frame(t1kk[,i])
  colnames(kkcumlist[[i-1]]) <-"status" 
  noinf <- filter(kkcumlist[[i-1]],status==0)
  noinf$rank <- 0
  kkcumlist[[i-1]] <- kkcumlist[[i-1]][kkcumlist[[i-1]]$status!=0,]
  kkcumlist[[i-1]] <- as.data.frame(kkcumlist[[i-1]])
  kkcumlist[[i-1]]$rank <- 1:nrow(kkcumlist[[i-1]])
  kkcumlist[[i-1]] <- bind_rows(noinf, kkcumlist[[i-1]])
  kkcumlist[[i-1]]$cumtot <- (kkcumlist[[i-1]]$rank/210)*100
}
kkcumarray <- array(as.numeric(unlist(kkcumlist)), dim=c(210,4,500))
kkcumquants <- as.data.frame(matrix(nrow=210, ncol = 3))
#take the quantiles across the sampled 500
for(j in 1:210){
  kkcumquants[j,] <- quantile(kkcumarray[j,4,1:500], probs=c(0.025, 0.5, 0.975))
}

ccacumquants[max(which(kkcumquants$V2==0.0000000)),]
ccacumquants[which(kkcumquants$V2==50.0000000),]
ccacumquants[which(kkcumquants$V2==10.0000000),]
gscorecumquants[max(which(kkcumquants$V2==0.0000000)),]
gscorecumquants[which(kkcumquants$V2==50.0000000),]
gscorecumquants[which(kkcumquants$V2==10.0000000),]

dev.new(height = 6, width = 12)
par(mfrow=c(1,2))
par(mar=c(6,6,4,1),mgp=c(3,1,0))

plot(NA,NA, xlim = c(0,100), ylim = c(0,100), xlab="cumulative % kato-katz infections", ylab="cumulative % cca+ infections", cex.axis=1)
lines(ccacumquants$V2~kkcumquants$V2, lwd = 2.5, col="#7570b3")
lines(ccacumquants$V1~kkcumquants$V1, lwd = 2.5, col="grey")
lines(ccacumquants$V3~kkcumquants$V3, lwd = 2.5, col="grey")
abline(v=50, col="lightgrey", lty=2)
abline(h=68.09524, col="lightgrey", lty=2)
abline(v=10, col="lightgrey", lty=2)
abline(h=28.09524, col="lightgrey", lty=2)

plot(NA,NA, xlim = c(0,100), ylim = c(0,100), xlab="cumulative % kato-katz infections", ylab="cumulative % G-Score infections", cex.axis=1)
lines(gscorecumquants$V2~kkcumquants$V2, lwd = 2.5, col="#d95f02")
lines(gscorecumquants$V1~kkcumquants$V1, lwd = 2.5, col="grey")
lines(gscorecumquants$V3~kkcumquants$V3, lwd = 2.5, col="grey")
abline(v=50, col="lightgrey", lty=2)
abline(h=61.90476, col="lightgrey", lty=2)
abline(v=10, col="lightgrey", lty=2)
abline(h=21.90476, col="lightgrey", lty=2)

dev.copy(pdf, "cumulativeinfections.pdf", height = 6, width = 12)
dev.off()
graphics.off()

#### Simulations ####
save(kkcca.model, file="kkccamodel.Rdata")
prevCCA <- as.numeric(combine.mcmc(kkcca.model$mcmc[,10:13])) 
shapes <- as.numeric(combine.mcmc(kkcca.model$mcmc[,2:5]))
rates <- as.numeric(combine.mcmc(kkcca.model$mcmc[,6:9]))
samples <- data.frame("prevalence" = prevCCA,"shapes" = shapes, "rates" = rates)
samples$means <- (samples$shapes/samples$rates)
samples$time <- c(rep("baseline",40000),rep("3 weeks",40000),rep("9 weeks",40000),rep("6 months",40000)) 

plot(NA,xlim = c(0,250),ylim = c(0,1), xlab = "", ylab = "",pch=19)
points((samples$shapes/samples$rates)*24,samples$prevalence,pch=19)

pdf("MeanvsPrevPoly.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
plot(NA,xlim = c(0,15),ylim = c(0,1),axes = F, xlab = "", ylab = "",pch=19)
#samples$means,samples$zeroinf into de NA for points
ids <- chull(samples$prevalence[as.factor(samples$time)=="baseline"],samples$means[as.factor(samples$time)=="baseline"])
polygon(samples$means[ids],samples$prevalence[ids],
        col=adjustcolor("darkgreen", alpha=.7), border = NA)
ids <- chull(samples$prevalence[as.factor(samples$time)=="3 weeks"],samples$means[as.factor(samples$time)=="3 weeks"])
polygon(samples$means[as.factor(samples$time)=="3 weeks"][ids],samples$prevalence[as.factor(samples$time)=="3 weeks"][ids],
        col=adjustcolor("black", alpha=.7), border = NA)
ids <- chull(samples$prevalence[as.factor(samples$time)=="9 weeks"],samples$means[as.factor(samples$time)=="9 weeks"])
polygon(samples$means[as.factor(samples$time)=="9 weeks"][ids],samples$prevalence[as.factor(samples$time)=="9 weeks"][ids],
        col=adjustcolor("darkblue", alpha=.7), border = NA)
ids <- chull(samples$prevalence[as.factor(samples$time)=="6 months"],samples$means[as.factor(samples$time)=="6 months"])
polygon(samples$means[as.factor(samples$time)=="6 months"][ids],samples$prevalence[as.factor(samples$time)=="6 months"][ids],
        col=adjustcolor("pink", alpha=.7), border = NA)
mtext("True prevalence",side=2,cex=1,line=1.2)
mtext("Average infection level",side=1,cex=1,line=1.2)

axis(1)
axis(2)

legend("bottomright", unique(samples$time), 
       col=c(adjustcolor("darkgreen", alpha=.7),adjustcolor("black", alpha=.7),adjustcolor("darkblue", alpha=.7)), 
       lty=c(1,1), bty='n',cex=0.75)

dev.off()

get.col <- function(alpha,col){
  return(adjustcolor(col,alpha))
}

do.plot <- function(timeid,col){
  z <- kde2d(samples$means[as.factor(samples$time)==timeid],samples$prevalence[as.factor(samples$time)==timeid], n=50)
  my.cols<-sapply(seq(0.2,1,by=.05),get.col,col)
  contour(z, drawlabels=FALSE, nlevels=length(my.cols), col=my.cols, add=TRUE)
}

pdf("MeanvsPrevCont.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
plot(NA,xlim = c(0,15),ylim = c(0,1),axes = F, xlab = "", ylab = "",pch=19)
#samples$means,samples$zeroinf into de NA for points
do.plot("baseline","darkgreen")
do.plot("3 weeks","black")
do.plot("9 weeks","darkblue")
do.plot("6 months","pink")

mtext("True prevalence",side=2,cex=1,line=1.2)
mtext("Average intensity of infection (epg)",side=1,cex=1,line=1.2)

xax <- c(0,5,10,15)
axis(1,at=xax,labels = xax*24)
axis(2)

legend("bottomright", unique(samples$time), 
       col=c(adjustcolor("darkgreen", alpha=.7),adjustcolor("black", alpha=.7),adjustcolor("darkblue", alpha=.7), adjustcolor("pink", alpha=.7)), 
       lty=c(1,1), bty='n',cex=0.75)

dev.off()



QuantilesSamples <- as.data.frame(apply(as.matrix(kkcca.model$mcmc[1]),2,quantile, probs=c(0.025,.5,0.95)))


# this doesn't do what I thought it would do - I thought the shape/rate would give me the mean of the egg counts
# but that can't be right because 
values <- data.frame("baselineprev"=QuantilesSamples$`prev[1]`, "threewkpev"=QuantilesSamples$`prev[2]`, 
                     "ninewkpev"=QuantilesSamples$`prev[3]`, "sixmnthpev"=QuantilesSamples$`prev[4]`, 
                     "blint"=QuantilesSamples$`tkksh[1]`/QuantilesSamples$`tkkrt[1]`, "thrint"=QuantilesSamples$`tkksh[2]`/QuantilesSamples$`tkkrt[2]`, 
                     "nineint"=QuantilesSamples$`tkksh[3]`/QuantilesSamples$`tkkrt[3]`, "sixint"=QuantilesSamples$`tkksh[4]`/QuantilesSamples$`tkkrt[4]`)
        
values$quants <- c("2.5%", "50%", "97.5%")    


sampleprev <- function(x){
  return(runif(1,min=lowerbound(x),max=upperbound(x)))
}

nDraws <- 1000
rho <- cor(samples$means,samples$rates)

library(mvtnorm)
Y <- rmvnorm(nDraws, sigma=matrix(c(1, rho, rho, 1), nrow=2)) # why is this a multivariate draw?
X <- pnorm(Y)
X[,1] <- X[,1]*(max(samples$means)-min(samples$means)) + min(samples$means)# this is rescaling 
X[,2] <- X[,2]*(max(samples$rates)-min(samples$rates)) + min(samples$rates)# this is rescaling

means <- X[,1]
rates <- X[,2]
shapes <- means * rates
prevalence <-  sapply( means ,sampleprev) #limits prev based on the mean

nRep = 1
nFECrep = 2
nInd = 1000

runRepeats <- function(index){
  return(apply(replicate(n=nRep,modelData(prevCCA=prevCCA[index],sh=shapes[index],rt=rates[index],nFECrep = nFECrep)),2,as.list))
}


modelData <- function(nInd=1000,prevCCA,sh,rt,nFECrep=2){
  
  k <- klogKKcca[runif(1,1,length(klogKKcca))]
  inter <- interceptKKcca[runif(1,1,length(interceptKKcca))]
  
  Status <- rbinom(nInd, size = 1, prob = prevCCA)
  Infectlvl <- rgamma(length(which(Status==1)),shape = sh, rate = rt)
  
  FECall <- matrix(rep(0,nInd*nFECrep),nrow=nFECrep)
  FECall[,which(Status==1)] <- sapply(Infectlvl, function(x) {return(rpois(n=nFECrep,lambda = x))})
  
  CCAall <- rep(0,nInd)
  CCAall[which(Status==0)] <- round(rtruncnorm(n=1, mean=0, a=0, b=1, sd=sqrt(1/3.093451)),0) # make truncated normal and change the 1 
  #logistic <- 4 / (1 + exp(-k*(Infectlvl-inter)))
  #CCAall[which(Status==1)] <- round(rtruncnorm(1, mean=0, a=0, sd=sqrt(1/3.093451)),0)
  #CCAall[which(Status==1)] <- sapply((4 / (1 + exp(-k*(Infectlvl-inter)))), function(x){round(rnorm(1, mean=x, 1/3.093451),0)})
  CCAall[which(Status==1)] <- sapply((4 / (1 + exp(-k*(Infectlvl-inter)))), function(x){round(rtruncnorm(n=1, mean=x, b=4, sd=sqrt(1/3.093451)),0)})
  
  #sapply(1:20, function(x){….})
  
  return(list(Status,FECall,CCAall,Infectlvl))
}

SimPopCCA <- t(sapply(1:length(shapes),runRepeats))

getTruePrevme <- function(exp){
  return(length(which(exp[[1]]==1))/nInd)
}

getFECPrevme <- function(exp){
  return(length(which(colMeans(matrix(as.numeric(exp[[2]]),ncol=nInd))!=0))/nInd)
}

getCCAPrevtpme <- function(exp){
  return(length(which(exp[[3]]!=0))/nInd)
}

getCCAPrevtnme <- function(exp){
  return(length(which(exp[[3]]>1))/nInd)
}

getMeanInflvlme <- function(exp){
  return(mean(c(exp[[4]],rep(0,length(exp[[1]]==0)))))
}

getMeanInflvlNozme <- function(exp){
  return(mean(c(exp[[4]])))
}

TruePrev <- matrix(as.numeric(lapply(SimPopCCA,getTruePrevme)),ncol=nRep)
TruePrevMean <- rowMeans(TruePrev)

FECPrev <- matrix(as.numeric(lapply(SimPopCCA,getFECPrevme)),ncol=nRep)
FECPrevMean <- rowMeans(FECPrev)

# not sure these are relevant because of the way we have done it 
CCAPrevtp <- matrix(as.numeric(lapply(SimPopCCA,getCCAPrevtpme)),ncol=nRep)
CCAPrevtpMean <- rowMeans(CCAPrevtp)

CCAPrevtn <- matrix(as.numeric(lapply(SimPopCCA,getCCAPrevtnme)),ncol=nRep)
CCAPrevtnMean <- rowMeans(CCAPrevtn)
####

MeansInflvl <- matrix(as.numeric(lapply(SimPopCCA,getMeanInflvlme)),ncol=nRep)
MeansInflvlmean <- rowMeans(MeansInflvl)

MeansInflvlNoz <- matrix(as.numeric(lapply(SimPopCCA,getMeanInflvlNozme)),ncol=nRep)
MeansInflvlNozmean <- rowMeans(MeansInflvlNoz)




