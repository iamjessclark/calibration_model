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
status <- model.output[,3:842]
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

prevKKcca <- density(c(kkcca.model$mcmc[[1]][,"prev"], kkcca.model$mcmc[[2]][,"prev"]))
rtnbKKcca <- density(c(kkcca.model$mcmc[[1]][,"rtnb"], kkcca.model$mcmc[[2]][,"rtnb"]))
klogKKcca <- density(c(kkcca.model$mcmc[[1]][,"k"], kkcca.model$mcmc[[2]][,"k"]))
interceptKKcca <- density(c(kkcca.model$mcmc[[1]][,"intercept"], kkcca.model$mcmc[[2]][,"intercept"]))

model.outputkkcca <- as.data.frame(as.matrix(as.mcmc(kkcca.model)))

# probability of infection at each time point #
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

# calibration of kk counts to cca # 
kk.kkccamod <- model.outputkkcca[845:5884]
kkccaT1C1 <- kk.kkccamod[,1:210]
kkccaT2C1 <- kk.kkccamod[,211:420]
kkccaT3C1 <- kk.kkccamod[,421:630]
kkccaT4C1 <- kk.kkccamod[,631:840]

kkccaT1C2 <- kk.kkccamod[,841:1050]
kkccaT2C2 <- kk.kkccamod[,1051:1260]
kkccaT3C2 <- kk.kkccamod[,1261:1470]
kkccaT4C2 <- kk.kkccamod[,1471:1680]

kkccaT1C3 <- kk.kkccamod[,1681:1890]
kkccaT2C3 <- kk.kkccamod[,1891:2100]
kkccaT3C3 <- kk.kkccamod[,2101:2310]
kkccaT4C3 <- kk.kkccamod[,2311:2520]

kkccaT1C4 <- kk.kkccamod[,2521:2730]
kkccaT2C4 <- kk.kkccamod[,2731:2940]
kkccaT3C4 <- kk.kkccamod[,2941:3150]
kkccaT4C4 <- kk.kkccamod[,3151:3360]

kkccaT1C5 <- kk.kkccamod[,3361:3570]
kkccaT2C5 <- kk.kkccamod[,3571:3780]
kkccaT3C5 <- kk.kkccamod[,3781:3990]
kkccaT4C5 <- kk.kkccamod[,3991:4200]

kkccaT1C6 <- kk.kkccamod[,4201:4410]
kkccaT2C6 <- kk.kkccamod[,4411:4620]
kkccaT3C6 <- kk.kkccamod[,4621:4830]
kkccaT4C6 <- kk.kkccamod[,4831:5040]

t1array <- list(kkccaT1C1,kkccaT1C2, kkccaT1C3,kkccaT1C4, kkccaT1C5, kkccaT1C6)
t1counts <- time.pointmeans(t2array)
rm(t1array)

t2array <- list(kkccaT2C1,kkccaT2C2, kkccaT2C3,kkccaT2C4, kkccaT2C5, kkccaT2C6)
t2counts <- time.pointmeans(t2array)
rm(t2array)

t3array <- list(kkccaT3C1,kkccaT3C2, kkccaT3C3,kkccaT3C4, kkccaT3C5, kkccaT3C6)
t3counts <- time.pointmeans(t3array)
rm(t3array)

t4array <- list(kkccaT4C1,kkccaT4C2, kkccaT4C3,kkccaT4C4, kkccaT4C5, kkccaT4C6)
t4counts <- time.pointmeans(t4array)
rm(t4array)

#### KK and G-Score Model ####
# remove the kkcca.model and model.outputkkcca objects as you don't have enough memory
# writing a csv doesn't work the model outputs will have to be unloaded before doing the next model
rm(kkcca.model)
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

#### Calibration Figures ####

