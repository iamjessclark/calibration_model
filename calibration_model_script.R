# Calibrating schistosomiasis diagnostics
# jessica clark 
# June 2020 

# packages #
source("diagnositic.calibration.model.R")
require(tidyverse)
require(runjags)
require(coda)

# load data #

kk.data <- read.csv("KK.4tp.csv") 
cca.data <- read.csv("cca.clean.csv")

# KK data for the model 
kk <- loadKKdata("KK.4tp.csv")
KKIDs <- getKKChildIDs("KK.4tp.csv")

# cca data for the model 
cca <- loadCCAdata("cca.clean.csv")
CCAIDs <- getCCAChildIDs("cca.clean.csv")

CID <- intersect(KKIDs, CCAIDs)
kk <- kk[match(CID, KKIDs),,]
cca <- cca[match(CID, CCAIDs),]

# g score data for the model 

ccagscore <- loadgscoredata("cca.clean.csv")
gscoreCCAIDS <- getCCAChildIDsgscore("cca.clean.csv") 

gscoreCID <- intersect(KKIDs,gscoreCCAIDS) 
ccagscore <- ccagscore[match(CID,gscoreCID),]

#### just KK model ####

katokatz <- CalModKK(dim(kk)[1], dim(kk)[2], dim(kk)[3], kk)
prevKK <- density(c(katokatz$mcmc[[1]][,"prev"], katokatz$mcmc[[2]][,"prev"]))
model.output <- as.data.frame(as.matrix(as.mcmc(katokatz)))
status <- model.output[,5:ncol(model.output)]
status.time.KK <- time.steps(status)

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

#### KK and CCA Model ####

kkcca.model <- CalModKKCCA(dim(kk)[1],dim(kk)[2],dim(kk)[3],kk,cca)
# monitoring prev, rtnb, k, intercept, Status, CCA, tKK
prevKKcca <- density(c(kkcca.model$mcmc[[1]][,"prev"], kkcca.model$mcmc[[2]][,"prev"]))
model.outputkkcca <- as.data.frame(as.matrix(as.mcmc(prevKKcca)))
status.kkcca <- model.outputkkcca[,5:844]
status.time.KKcca <- time.steps(status.kkcca)

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

#### KK and G-Score Model ####

kkGscore.model <- CalModKKGScore(dim(kk)[1],dim(kk)[2],dim(kk)[3],kk,ccagscore)

prevKKGScore <- density(c(kkGscore.model$mcmc[[1]][,"prev"], kkGscore.model$mcmc[[2]][,"prev"]))
outputkkGscore <- as.data.frame(as.matrix(as.mcmc(kkGscore.model)))
status.kkgscore <- outputkkGscore[,5:844]
status.time.KKgscore <- time.steps(status.kkgscore)

pdf("probinfKKgscore.pdf",width=3,height=8)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfrow=c(1,1))

prob.inf.kk <- qplot(value, data = status.time.KKgscore, geom = "histogram", fill=time, binwidth=0.02) +
  scale_fill_manual(values = c("#e7298a", "#1b9e77", "#d95f02", "#7570b3"))+
  geom_text(data=anno, aes(x = xcoord, y = ycoord, label = label))+
  facet_grid(time~.)+
  theme_bw()+coord_cartesian(ylim=c(0, 200), xlim=c(0,1))+theme(legend.position = "none")+
  ylab("Number of Hosts") + xlab("Probability of Being Infected")

prob.inf.kk
dev.off()