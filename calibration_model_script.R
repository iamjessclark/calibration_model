# Calibrating schistosomiasis diagnostics
# jessica clark 
# June 2020 

# packages #
source("diagnositic.calibration.model.R")
require(data.table)
require(reshape2)
require(tidyverse)
require(runjags)
require(coda)
require(viridis)

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

# get your random sample to make the confidence intervals from 
kintkkcca <- as.data.frame(cbind(c(kkcca.model$mcmc[[1]][,"k"], kkcca.model$mcmc[[2]][,"k"]),c(kkcca.model$mcmc[[1]][,"intercept"], kkcca.model$mcmc[[2]][,"intercept"])))
kkcca.sample <- kintkkcca[sample.int(nrow(kintkkcca), 100, replace=FALSE, prob=NULL),]

# make model output into a df so that you can get the things you want out of it!
model.outputkkcca <- as.data.frame(as.matrix(as.mcmc(kkcca.model)))
status.kkcca <- model.outputkkcca[,5:844]
status.time.KKcca <- time.steps(status.kkcca)

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

#### KK and G-Score Model ####

# remove the kkcca.model and model.outputkkcca objects as you don't have enough memory
# writing a csv doesn't work the model outputs will have to be unloaded before doing the next model
rm(kkcca.model)
kkGscore.model <- CalModKKGScore(dim(kk)[1],dim(kk)[2],dim(kk)[3],kk,ccagscore)

prevKKGScore <- density(c(kkGscore.model$mcmc[[1]][,"prev"], kkGscore.model$mcmc[[2]][,"prev"]))
klogKKgscore <- density(c(kkGscore.model$mcmc[[1]][,"k"], kkGscore.model$mcmc[[2]][,"k"]))
interceptKKgscore <- density(c(kkGscore.model$mcmc[[1]][,"intercept"], kkGscore.model$mcmc[[2]][,"intercept"]))

#random samples to make CI from
kintkkgscore <- as.data.frame(cbind(c(kkGscore.model$mcmc[[1]][,"k"], kkGscore.model$mcmc[[2]][,"k"]),c(kkGscore.model$mcmc[[1]][,"intercept"], kkGscore.model$mcmc[[2]][,"intercept"])))
kkgscore.sample <- kintkkgscore[sample.int(nrow(kintkkgscore), 100, replace=FALSE, prob=NULL),]

# make the model into a df for plotting
outputkkGscore <- as.data.frame(as.matrix(as.mcmc(kkGscore.model)))
status.kkgscore <- outputkkGscore[,5:844]
status.time.KKgscore <- time.steps(status.kkgscore)

# logistic curve figure 

logcurvelistkkcca <- list()
logcurvelistkkgscore <- list()
x=seq(from=1, to=100, by=1)

for(i in 1:100){
  logcurvelistkkcca[[i]] <- as.data.frame(logcurve(x,kintkkcca[i,1],kintkkcca[i,2]))
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

for(i in 1:100){
  logcurvelistkkgscore[[i]] <- as.data.frame(logcurve(x,kintkkgscore[i,1],kintkkgscore[i,2]))
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

plot(NA,NA, xlim = c(0,100), ylim = c(0,1), xlab="Kato-Katz Count", ylab="f(x)", cex.axis=1)
lines(kkccaquantiles$`50%`, col="orange", lwd=2, lty = 1)
polygon(x=polygonx, y=polygony, col=adjustcolor("grey", alpha.f=0.4), border=NA)
lines(kkgscorequantiles$`2.5%`, col="darkred", lwd=1) 
polygon(x=polygonxgscore, y=polygonygscore, col=adjustcolor("grey", alpha.f=0.4), border=NA)

add_legend("top", legend=c(c("KK & CCA (normal)", "KK & G-Score")), lty=1, lwd=2,
           col=c("orange","darkred"),
           horiz=TRUE, bty='n', cex=1)
dev.copy(pdf, "logcurve.confint.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# intercept/k posterior plots #
dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
pdf("kposts.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(klogKKgscore$x,klogKKgscore$y/max(klogKKgscore$y),xlim=c(0,100),ylim=c(0,1),type = "l", ylab = "",xlab="", col="blue")
plot(klogKKcca$x,klogKKcca$y/max(klogKKcca$y),xlim=c(0,100),ylim=c(0,1),type = "l", ylab = "",xlab="", col="red")

mtext("Scaled density",side=2,cex=1,line=1.2)
mtext("logistic function intercept ",side=1,cex=1,line=1.2)
dev.off()
graphics.off()

pdf("probinfKKgscore.pdf",width=3,height=8)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfrow=c(1,1))

# probability of being infected at each time point 

prob.inf.kkgscore <- qplot(value, data = status.time.KKgscore, geom = "histogram", fill=time, binwidth=0.02) +
  scale_fill_manual(values = c("#e7298a", "#1b9e77", "#d95f02", "#7570b3"))+
  geom_text(data=anno, aes(x = xcoord, y = ycoord, label = label))+
  facet_grid(time~.)+
  theme_bw()+coord_cartesian(ylim=c(0, 200), xlim=c(0,1))+theme(legend.position = "none")+
  ylab("Number of Hosts") + xlab("Probability of Being Infected")

prob.inf.kk
dev.off()

#### raw data plots ####

# restructure for categorical
raw.cca.data <- cca.data
raw.cca.data$Poppy <- as.character(raw.cca.data$Poppy)
raw.cca.data$Poppy[which(raw.cca.data$Poppy=="T")] <- "Trace"
raw.cca.data$Poppy[which(raw.cca.data$Poppy==3)]<-"+++"
raw.cca.data$Poppy[which(raw.cca.data$Poppy==2)]<- "++"
raw.cca.data$Poppy[which(raw.cca.data$Poppy==1)]<- "+"
raw.cca.data$Poppy[which(raw.cca.data$Poppy==0.5)]<- "Trace"
raw.cca.data$Poppy[which(raw.cca.data$Poppy==0)]<- "Negative (0)"
raw.cca.data$Poppy <- factor(raw.cca.data$Poppy)
raw.cca.data <- raw.cca.data %>% filter(Poppy!="-")
raw.cca.data$Poppy <- factor(raw.cca.data$Poppy)
raw.cca.data <- raw.cca.data %>% filter(Poppy!="")
raw.cca.data$Poppy <- factor(raw.cca.data$Poppy, levels=c("Negative (0)", "Trace", "+", "++","+++"))
raw.cca.data$gscore <- as.integer(raw.cca.data$gscore)

raw.cca.data$dateN <- NA
raw.cca.data$dateN[which(raw.cca.data$date_on_tube=="25/09/2017" | raw.cca.data$date_on_tube=="26/09/2017" | raw.cca.data$date_on_tube=="27/09/2017"
                         | raw.cca.data$date_on_tube=="28/09/2017" | raw.cca.data$date_on_tube=="29/09/2017" | raw.cca.data$date_on_tube=="02/10/2017")] <- "Pre-T"

raw.cca.data$dateN[which(raw.cca.data$date_on_tube=="23/10/2017" | raw.cca.data$date_on_tube=="24/10/2017" | raw.cca.data$date_on_tube=="25/10/2017"
                         | raw.cca.data$date_on_tube=="26/10/2017" | raw.cca.data$date_on_tube=="27/10/2017")] <- "3 weeks"

raw.cca.data$dateN[which(raw.cca.data$date_on_tube=="04/12/2017" | raw.cca.data$date_on_tube=="05/12/2017")] <- "9 weeks"

raw.cca.data$dateN[which(raw.cca.data$date_on_tube=="01/03/2018" | raw.cca.data$date_on_tube=="05/03/2018" | raw.cca.data$date_on_tube=="06/03/2018"
                         | raw.cca.data$date_on_tube=="07/03/2018" | raw.cca.data$date_on_tube=="08/03/2018" | raw.cca.data$date_on_tube=="09/03/2018")] <- "6 months"

ccadates <- raw.cca.data[,c("CID", "Poppy", "gscore", "dateN")]
ccadates$dateN <- as.factor(ccadates$dateN)

# sort out KK data for plot 

kkxgscore.data <- kk.data
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

corr.diag.plot <- ggplot(data=mergekkcca, aes(x=gscore, y=Poppy))+
  geom_jitter(aes(colour=intcat), size=2)+
  facet_grid(dateN~.)+
  labs(colour="Intensity Categories", y="CCA Standard", x="G-Score")+
  scale_colour_viridis(discrete=TRUE) +
  ggsave("corr.diag.plot.raw.pdf")
  
#### Sensitivity and specificity ROC plot ####
mergekkcca$CID <- factor(mergekkcca$CID)
mergekkcca$gscore <- as.numeric(mergekkcca$gscore)
mergekkcca <- mergekkcca %>% mutate(kkpos = ifelse(mean.epg>0, 1, 0), 
                                    ccaTpos = ifelse(Poppy == "Trace" |Poppy ==  "+" | Poppy == "++" |Poppy == "+++", 1,0),
                                    cca1pos = ifelse(Poppy == "+" |Poppy ==  "++" |Poppy == "+++", 1,0),
                                    cca2pos = ifelse(Poppy == " ++" |Poppy == "+++", 1,0),
                                    cca3pos = ifelse(Poppy == "+++", 1,0),
                                    gscore2pos = ifelse(gscore>=2, 1, 1), 
                                    gscore3pos = ifelse(gscore>=3, 1, 0), 
                                    gscore4pos = ifelse(gscore>=4, 1, 0), 
                                    gscore5pos = ifelse(gscore>=5, 1, 0), 
                                    gscore6pos = ifelse(gscore>=6, 1, 0), 
                                    gscore7pos = ifelse(gscore>=7, 1, 0), 
                                    gscore8pos = ifelse(gscore>=8, 1, 0),
                                    gscore9pos = ifelse(gscore>=9, 1, 0), 
                                    gscore10pos = ifelse(gscore==10, 1, 0))


#kato katz model #
kkt1 <- as.data.frame(status[,1:210])
kkt2 <- as.data.frame(status[,211:420])
kkt3 <- as.data.frame(status[,421:630])
kkt4 <- as.data.frame(status[,631:840])

colnames(kkt1) <- CID
colnames(kkt2) <- CID
colnames(kkt3) <- CID
colnames(kkt4) <- CID

# t1sample 500 random particals 
kkt1status.sample <- kkt1[sample.int(nrow(kkt1), 500, replace=FALSE, prob=NULL),]
kkt1status.sample <- t(kkt1status.sample)
kkt1status.sample <- as.data.frame(cbind(rownames(kkt1status.sample), data.frame(kkt1status.sample, row.names=NULL)))
kkt1status.sample <- kkt1status.sample %>% rename(CID=`rownames(kkt1status.sample)`)

# join to the kkpos
mergepretkk <- mergekkcca %>% select(dateN, CID, kkpos) %>%
                filter(dateN=="Pre-T")%>%
                  right_join(kkt1status.sample, by="CID")

# remove duplicates
mergepretkk <- mergepretkk[-which(is.na(mergepretkk$kkpos)==T),]
mergepretkk <- mergepretkk[!duplicated(mergepretkk[,"CID"]),]

# get the sensitivity/ specificity 
for(i in 1:nrow(mergepretkk)){
  for(c in 4:ncol(mergepretkk)){
  if(mergepretkk[i,3]==1 & mergepretkk[i,c]==1){
    mergepretkk[i,c] <- "TP"
  } else if(mergepretkk[i,3]==1 & mergepretkk[i,c]==0){
    mergepretkk[i,c] <- "FP"
  } else if(mergepretkk[i,3]==0 & mergepretkk[i,c]==0){
    mergepretkk[i,c] <- "TN" 
  } else {
    mergepretkk[i,c] <- "FN" 
    }
  }
}

mergepretkk <- mergepretkk %>% mutate_if(is.character, as.factor)
pretkklist <- list()

# calculate how often each one occurs and get proportions 
for(i in 4:ncol(mergepretkk)){
  pretkklist[[i-3]] <-  mergepretkk %>% group_by(mergepretkk[,i]) %>% tally() %>%
                        mutate(proportion = n/sum(n))
}

test.specs.kk <- list()
for(i in 1:length(pretkklist)){
  test.specs.kk[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
  colnames(test.specs.kk[[i]]) <- c("TPR", "FPR")
  test.specs.kk[[i]][1,1] <- pretkklist[[i]][4,2]/(pretkklist[[i]][4,2]+pretkklist[[i]][1,2])
  test.specs.kk[[i]][1,2] <- pretkklist[[i]][2,2]/(pretkklist[[i]][3,2]+pretkklist[[i]][2,2])
}

roc.pretkk <- rbindlist(test.specs.kk)
roc.pretkk <- roc.pretkk %>% mutate(time = "pre-T")%>%mutate_if(is.character, as.factor)

# t2sample 500 random particals 
kkt2status.sample <- kkt2[sample.int(nrow(kkt2), 500, replace=FALSE, prob=NULL),]
kkt2status.sample <- t(kkt2status.sample)
kkt2status.sample <- as.data.frame(cbind(rownames(kkt2status.sample), data.frame(kkt2status.sample, row.names=NULL)))
kkt2status.sample <- kkt2status.sample %>% rename(CID=`rownames(kkt2status.sample)`)

mergekk3wk <- mergekkcca %>% select(dateN, CID, kkpos) %>%
  filter(dateN=="3 weeks")%>%
  right_join(kkt2status.sample, by="CID")
mergekk3wk <- mergekk3wk[-which(is.na(mergekk3wk$kkpos)==T),]
mergekk3wk <- mergekk3wk[!duplicated(mergekk3wk[,"CID"]),]

# get the sensitivity/ specificity 
for(i in 1:nrow(mergekk3wk)){
  for(c in 4:ncol(mergekk3wk)){
    if(mergekk3wk[i,3]==1 & mergekk3wk[i,c]==1){
      mergekk3wk[i,c] <- "TP"
    } else if(mergekk3wk[i,3]==1 & mergekk3wk[i,c]==0){
      mergekk3wk[i,c] <- "FP"
    } else if(mergekk3wk[i,3]==0 & mergekk3wk[i,c]==0){
      mergekk3wk[i,c] <- "TN" 
    } else {
      mergekk3wk[i,c] <- "FN" 
    }
  }
}

mergekk3wk <- mergekk3wk %>% mutate_if(is.character, as.factor)
wk3kklist <- list()

# calculate how often each one occurs and get proportions 
for(i in 4:ncol(mergekk3wk)){
  wk3kklist[[i-3]] <-  mergekk3wk %>% group_by(mergekk3wk[,i]) %>% tally() %>%
    mutate(proportion = n/sum(n))
}

test.specs.wk3kk <- list()
for(i in 1:length(wk3kklist)){
  test.specs.wk3kk[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
  colnames(test.specs.wk3kk[[i]]) <- c("TPR", "FPR")
  test.specs.wk3kk[[i]][1,1] <- wk3kklist[[i]][4,2]/(wk3kklist[[i]][4,2]+wk3kklist[[i]][1,2])
  test.specs.wk3kk[[i]][1,2] <- wk3kklist[[i]][2,2]/(wk3kklist[[i]][3,2]+wk3kklist[[i]][2,2])
}

roc.3wk <- rbindlist(test.specs.wk3kk)
roc.3wk <- roc.3wk %>% mutate(time = "3 weeks")%>%mutate_if(is.character, as.factor)
# t3sample 500 random particals 
kkt3status.sample <- kkt3[sample.int(nrow(kkt3), 500, replace=FALSE, prob=NULL),]
kkt3status.sample <- t(kkt3status.sample)
kkt3status.sample <- as.data.frame(cbind(rownames(kkt3status.sample), data.frame(kkt3status.sample, row.names=NULL)))
kkt3status.sample <- kkt3status.sample %>% rename(CID=`rownames(kkt3status.sample)`)

mergekk9wk <- mergekkcca %>% select(dateN, CID, kkpos) %>%
  filter(dateN=="3 weeks")%>%
  right_join(kkt3status.sample, by="CID")

mergekk9wk <- mergekk9wk[-which(is.na(mergekk9wk$kkpos)==T),]
mergekk9wk <- mergekk9wk[!duplicated(mergekk9wk[,"CID"]),]

# get the sensitivity/ specificity 
for(i in 1:nrow(mergekk9wk)){
  for(c in 4:ncol(mergekk9wk)){
    if(mergekk9wk[i,3]==1 & mergekk9wk[i,c]==1){
      mergekk9wk[i,c] <- "TP"
    } else if(mergekk9wk[i,3]==1 & mergekk9wk[i,c]==0){
      mergekk9wk[i,c] <- "FP"
    } else if(mergekk9wk[i,3]==0 & mergekk9wk[i,c]==0){
      mergekk9wk[i,c] <- "TN" 
    } else {
      mergekk9wk[i,c] <- "FN" 
    }
  }
}

mergekk9wk <- mergekk9wk %>% mutate_if(is.character, as.factor)
wk9kklist <- list()

# calculate how often each one occurs and get proportions 
for(i in 4:ncol(mergekk9wk)){
  wk9kklist[[i-3]] <-  mergekk9wk %>% group_by(mergekk9wk[,i]) %>% tally() %>%
    mutate(proportion = n/sum(n))
}

test.specs.wk9kk <- list()
for(i in 1:length(wk9kklist)){
  test.specs.wk9kk[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
  colnames(test.specs.wk9kk[[i]]) <- c("TPR", "FPR")
  test.specs.wk9kk[[i]][1,1] <- wk9kklist[[i]][4,2]/(wk9kklist[[i]][4,2]+wk9kklist[[i]][1,2])
  test.specs.wk9kk[[i]][1,2] <- wk9kklist[[i]][2,2]/(wk9kklist[[i]][3,2]+wk9kklist[[i]][2,2])
}

roc.9wk <- rbindlist(test.specs.wk9kk)
roc.9wk <- roc.9wk %>% mutate(time = "9 weeks")%>%mutate_if(is.character, as.factor)

# t3sample 500 random particals 
kkt4status.sample <- kkt4[sample.int(nrow(kkt4), 500, replace=FALSE, prob=NULL),]
kkt4status.sample <- t(kkt4status.sample)
kkt4status.sample <- as.data.frame(cbind(rownames(kkt4status.sample), data.frame(kkt4status.sample, row.names=NULL)))
kkt4status.sample <- kkt4status.sample %>% rename(CID=`rownames(kkt4status.sample)`)

mergekk6m <- mergekkcca %>% select(dateN, CID, kkpos) %>%
  filter(dateN=="6 months")%>%
  right_join(kkt4status.sample, by="CID")

# remove duplicates
mergekk6m <- mergekk6m[-which(is.na(mergekk6m$kkpos)==T),]
mergekk6m <- mergekk6m[!duplicated(mergekk6m[,"CID"]),]

# get the sensitivity/ specificity 
for(i in 1:nrow(mergekk6m)){
  for(c in 4:ncol(mergekk6m)){
    if(mergekk6m[i,3]==1 & mergekk6m[i,c]==1){
      mergekk6m[i,c] <- "TP"
    } else if(mergekk6m[i,3]==1 & mergekk6m[i,c]==0){
      mergekk6m[i,c] <- "FP"
    } else if(mergekk6m[i,3]==0 & mergekk6m[i,c]==0){
      mergekk6m[i,c] <- "TN" 
    } else {
      mergekk6m[i,c] <- "FN" 
    }
  }
}

mergekk6m <- mergekk6m %>% mutate_if(is.character, as.factor)
sixmkklist <- list()

# calculate how often each one occurs and get proportions 
for(i in 4:ncol(mergekk6m)){
  sixmkklist[[i-3]] <-  mergekk6m %>% group_by(mergekk6m[,i]) %>% tally() %>%
    mutate(proportion = n/sum(n))
}

test.specs.sixmkk <- list()
for(i in 1:length(sixmkklist)){
  test.specs.sixmkk[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
  colnames(test.specs.sixmkk[[i]]) <- c("TPR", "FPR")
  test.specs.sixmkk[[i]][1,1] <- sixmkklist[[i]][4,2]/(sixmkklist[[i]][4,2]+sixmkklist[[i]][1,2])
  test.specs.sixmkk[[i]][1,2] <- sixmkklist[[i]][2,2]/(sixmkklist[[i]][3,2]+sixmkklist[[i]][2,2])
}

roc.6m <- rbindlist(test.specs.sixmkk)
roc.6m <- roc.6m %>% mutate(time = "6 months")%>%mutate_if(is.character, as.factor)

kk.spec <- bind_rows(roc.pretkk, roc.3wk, roc.9wk, roc.6m)%>% mutate_if(is.character, as.factor)
kk.spec$time <- factor(kk.spec$time, levels=c("pre-T", "3 weeks", "9 weeks", "6 months"))

kk.spec.qtiles <- kk.spec %>% group_by(time) %>%
                  summarise_at(vars(TPR, FPR),funs(!!!p_funs))

kk.spec.plot <- ggplot()+
  geom_jitter(kk.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  ylab("True Positive Rate") + xlab("False Positive Rate")+
  ggsave("kkspecplot.pdf")

#### CCA Model T positive ####

condition.data.frame <- mergekkcca %>% select(dateN, CID, ccaTpos)
ccaTpos.spec <- diag_spec(status.kkcca, condition.data.frame) 
ccaTpos.spec$time <- factor(ccaTpos.spec$time, levels=c("pre-T", "3 weeks", "9 weeks", "6 months"))
specsensccaT <- ccaTpos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="trace")

condition.data.frame <- mergekkcca %>% select(dateN, CID, cca1pos)
cca1pos.spec <- diag_spec(status.kkcca, condition.data.frame) 
specsenscca1 <- cca1pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs))%>% 
  mutate(threshold="+")

condition.data.frame <- mergekkcca %>% select(dateN, CID, cca2pos)
cca2pos.spec <- diag_spec(status.kkcca, condition.data.frame) 
specsenscca2 <- cca2pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs))%>% 
  mutate(threshold="++")

condition.data.frame <- mergekkcca %>% select(dateN, CID, cca3pos)
cca3pos.spec <- diag_spec(status.kkcca, condition.data.frame) 
specsenscca3 <- cca3pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs))%>% 
  mutate(threshold="+++")

ccasensspecquant <- bind_rows(specsensccaT,specsenscca1, specsenscca2, specsenscca2 )%>%
  mutate_if(is.character, as.factor)


CCAT.spec.plot <- ggplot()+
  geom_point(ccaTpos.spec, mapping=aes(x=FPR , y=TPR))+
  geom_point(cca1pos.spec, mapping=aes(x=FPR , y=TPR))+
  geom_point(cca2pos.spec, mapping=aes(x=FPR , y=TPR))+
  geom_point(cca3pos.spec, mapping=aes(x=FPR , y=TPR))+
  ylab("True Positive Rate") + xlab("False Positive Rate")


CCAT.spec.plot <- ggplot()+
  geom_point(ccasensspecquant, mapping=aes(x=`FPR_50%` , y=`TPR_50%`, colour=threshold, shape=time))+
  geom_line()+
  ylab("True Positive Rate") + xlab("False Positive Rate")

polygonx <- c(ccasensspecquant$`FPR_50%`, rev(ccasensspecquant$`FPR_50%`))
polygony <- c(ccasensspecquant$`TPR_2.5%`, rev(ccasensspecquant$`TPR_97.5%`))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="False Positive Rate", ylab="True Positive Rate", cex.axis=1)
lines(ccasensspecquant$`FPR_50%`~ccasensspecquant$`TPR_50%`, col="orange", lwd=2, lty = 1)
polygon(x=polygonx, y=polygony, col=adjustcolor("grey", alpha.f=0.4), border=NA)

add_legend("top", legend=c(c("KK & CCA (normal)", "KK & G-Score")), lty=1, lwd=2,
           col=c("orange","darkred"),
           horiz=TRUE, bty='n', cex=1)
dev.copy(pdf, "logcurve.confint.pdf", height = 6, width = 6)
dev.off()
graphics.off()


CCAT.spec.plot <- ggplot()+
  geom_jitter(gscore2pos.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  geom_jitter(gscore3pos.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  geom_jitter(gscore4pos.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  geom_jitter(gscore5pos.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  geom_jitter(gscore6pos.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  geom_jitter(gscore7pos.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  geom_jitter(gscore8pos.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  geom_jitter(gscore9pos.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  geom_jitter(gscore10pos.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  ylab("True Positive Rate") + xlab("False Positive Rate")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore2pos)
gscore2pos.spec <- diag_spec(status.kkgscore, condition.data.frame)
specsensgscore2 <- gscore2pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score2")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore3pos)
gscore3pos.spec <- diag_spec(status.kkgscore, condition.data.frame) 
specsensgscore3 <- gscore3pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score3")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore4pos)
gscore4pos.spec <- diag_spec(status.kkgscore, condition.data.frame)
specsensgscore4 <- gscore4pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score4")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore5pos)
gscore5pos.spec <- diag_spec(status.kkgscore, condition.data.frame)
specsensgscore5 <- gscore5pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score5")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore6pos)
gscore6pos.spec <- diag_spec(status.kkgscore, condition.data.frame) 
specsensgscore6 <- gscore6pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score6")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore7pos)
gscore7pos.spec <- diag_spec(status.kkgscore, condition.data.frame)
specsensgscore7 <- gscore7pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score7")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore8pos)
gscore8pos.spec <- diag_spec(status.kkgscore, condition.data.frame) 
specsensgscore8 <- gscore8pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score8")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore9pos)
gscore9pos.spec <- diag_spec(status.kkgscore, condition.data.frame)
specsensgscore9 <- gscore9pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score9")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore10pos)
gscore10pos.spec <- diag_spec(status.kkgscore, condition.data.frame) 
specsensgscore10 <- gscore10pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score10")

gscoresensspecquant <- bind_rows(specsensgscore2, specsensgscore3, specsensgscore4,
                                 specsensgscore5, specsensgscore6, specsensgscore7, specsensgscore8, 
                                 specsensgscore9, specsensgscore10)%>% 
                        mutate(time=factor(time, levels=c("pre-T", "3 weeks", "9 weeks", "6 months")))%>%
                        mutate(threshold=factor(threshold, levels=c("G-Score2","G-Score3","G-Score4","G-Score5","G-Score6","G-Score7","G-Score8","G-Score9","G-Score10")))

gscoresensspecall <- bind_rows(gscore2pos.spec, gscore3pos.spec,gscore4pos.spec, gscore5pos.spec, gscore6pos.spec, gscore7pos.spec, gscore8pos.spec, gscore9pos.spec, gscore10pos.spec)

gscore.spec.plot <- ggplot()+
  geom_point(gscoresensspecquant, mapping=aes(x=`FPR_50%` , y=`TPR_50%`, colour=threshold))+
  geom_line(data=gscoresensspecall, aes(x=FPR, y=TPR))+
  geom_ribbon(data=gscoresensspecall, aes(ymin=TPR, ymax=))+
  facet_wrap(time~.)+
  ylab("True Positive Rate") + xlab("False Positive Rate")

