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

add_legend("top", legend=c(c("Plus (+)", "G-Score")), lty=1, lwd=2,
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

prob.inf.kkgscore
dev.off()

status.time.KKcca <- status.time.KKcca %>% mutate(diagnostic="Plus(+)")
status.time.KKgscore <- status.time.KKgscore %>% mutate(diagnostic="G-Score")
status.time.KK <- status.time.KK %>% mutate(diagnostic="Kato-Katz")

probinfdf <- bind_rows(status.time.KKcca,status.time.KKgscore, status.time.KK )%>%
  mutate_if(is.character, as.factor)

prob.inf.all.plot <- ggplot(data=probinfdf)+
  geom_histogram(bins=10, aes(x=value, fill=diagnostic), position="dodge")+
  facet_grid(time~.)+
  scale_fill_viridis_d(option = "plasma")+
  ylab("Number of Children") + xlab("Probability of Being Infected")+
  theme_bw()+
  ggsave("prob.inf.all.pdf")

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
  labs(colour="Intensity Categories", y="Plus (+) Score", x="G-Score")+
  scale_colour_viridis(discrete=TRUE) +
  theme_bw()+
  scale_x_continuous(breaks=c(1:10,1))+
  ggsave("corr.diag.plot.raw.pdf")
  
#### Sensitivity and specificity ROC plot ####
mergekkcca$CID <- factor(mergekkcca$CID)
mergekkcca$gscore <- as.numeric(mergekkcca$gscore)
mergekkcca <- mergekkcca %>% mutate(kkpos = ifelse(mean.epg>0, 1, 0), 
                                    ccanegpos = 1,
                                    ccaTpos = ifelse(Poppy == "Trace" |Poppy ==  "+" | Poppy == "++" |Poppy == "+++", 1,0),
                                    cca1pos = ifelse(Poppy == "+" |Poppy ==  "++" |Poppy == "+++", 1,0),
                                    cca2pos = ifelse(Poppy == " ++" |Poppy == "+++", 1,0),
                                    cca3pos = ifelse(Poppy == "+++", 1,0),
                                    gscore1pos = ifelse(gscore>=1, 1, 0), 
                                    gscore2pos = ifelse(gscore>=2, 1, 0), 
                                    gscore3pos = ifelse(gscore>=3, 1, 0), 
                                    gscore4pos = ifelse(gscore>=4, 1, 0), 
                                    gscore5pos = ifelse(gscore>=5, 1, 0), 
                                    gscore6pos = ifelse(gscore>=6, 1, 0), 
                                    gscore7pos = ifelse(gscore>=7, 1, 0), 
                                    gscore8pos = ifelse(gscore>=8, 1, 0),
                                    gscore9pos = ifelse(gscore>=9, 1, 0), 
                                    gscore10pos = ifelse(gscore==10, 1, 0))


#kato katz model #

kk.sample <- status[sample.int(nrow(status), 500, replace=FALSE, prob=NULL),]
kkt1 <- as.data.frame(kk.sample[,1:210])
kkt2 <- as.data.frame(kk.sample[,211:420])
kkt3 <- as.data.frame(kk.sample[,421:630])
kkt4 <- as.data.frame(kk.sample[,631:840])

colnames(kkt1) <- CID
colnames(kkt2) <- CID
colnames(kkt3) <- CID
colnames(kkt4) <- CID

kkt1 <- t(kkt1)
kkt2 <- t(kkt2)
kkt3 <- t(kkt3)
kkt4 <- t(kkt4)

# this turns row names into the column of CIDs
kkt1 <- as.data.frame(cbind(rownames(kkt1), data.frame(kkt1, row.names=NULL)))
kkt1 <- kkt1 %>% rename(CID=`rownames(kkt1)`)

mergepretkk <- mergekkcca %>% select(dateN, CID, kkpos) %>%
  right_join(kkt1, by="CID")%>%
  filter(dateN=="Pre-T")
# remove duplicates

kkt2 <- as.data.frame(cbind(rownames(kkt2), data.frame(kkt2, row.names=NULL)))
kkt2 <- kkt2 %>% rename(CID=`rownames(kkt2)`)

mergekk3wk <- mergekkcca %>% select(dateN, CID, kkpos) %>%
  right_join(kkt2, by="CID")%>%
  filter(dateN=="3 weeks")
# remove duplicates

kkt3 <- as.data.frame(cbind(rownames(kkt3), data.frame(kkt3, row.names=NULL)))
kkt3 <- kkt3 %>% rename(CID=`rownames(kkt3)`)

mergekk9wk <- mergekkcca %>% select(dateN, CID, kkpos) %>%
  right_join(kkt3, by="CID")%>%
  filter(dateN=="9 weeks")
# remove duplicates

kkt4 <- as.data.frame(cbind(rownames(kkt4), data.frame(kkt4, row.names=NULL)))
kkt4 <- kkt4 %>% rename(CID=`rownames(kkt4)`)

mergekk6m <- mergekkcca %>% select(dateN, CID, kkpos) %>%
  right_join(kkt4, by="CID")%>%
  filter(dateN=="6 months")
# remove duplicates


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

table(factor(mergepretkk$X9353))

for(i in 4:ncol(mergepretkk)){
  pretkklist[[i-3]] <-  mergepretkk %>% group_by(mergepretkk[,i]) %>% tally() %>%
                        mutate(proportion = n/sum(n))
}

test.specs.t1 <- list()
for(i in 1:length(pretkklist)){
  test.specs.t1[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
  colnames(test.specs.t1[[i]]) <- c("TPR", "FPR")
  
  if("TP" %in% pretkklist[[i]]$ `mergepretkk[, i]` == TRUE && "FN" %in% pretkklist[[i]]$ `mergepretkk[, i]` == TRUE){
    test.specs.t1[[i]][1,1] <- pretkklist[[i]][[which(pretkklist[[i]]=="TP"),2]]/(pretkklist[[i]][[which(pretkklist[[i]]=="TP"),2]]+pretkklist[[i]][[which(pretkklist[[i]]=="FN"),2]])
  } else if("TP" %in% pretkklist[[i]]$ `mergepretkk[, i]`  == TRUE && "FN" %in% pretkklist[[i]]$ `mergepretkk[, i]` ==FALSE){
    test.specs.t1[[i]][1,1] <- pretkklist[[i]][[which(pretkklist[[i]]=="TP"),2]]/(pretkklist[[i]][[which(pretkklist[[i]]=="TP"),2]]+0)
  } else {
    test.specs.t1[[i]][1,1] <- NA
  }
  
  if("FP" %in% pretkklist[[i]]$ `mergepretkk[, i]` == TRUE && "TN" %in% pretkklist[[i]]$ `mergepretkk[, i]` ==TRUE){
    test.specs.t1[[i]][1,2] <- (pretkklist[[i]][[which(pretkklist[[i]]=="FP"),2]]/(pretkklist[[i]][[which(pretkklist[[i]]=="TN"),2]]+pretkklist[[i]][[which(pretkklist[[i]]=="FP"),2]]))
  } else if("FP" %in% pretkklist[[i]]$ `mergepretkk[, i]`  == TRUE && "TN" %in% pretkklist[[i]]$ `mergepretkk[, i]` ==FALSE){
    test.specs.t1[[i]][1,2] <- (pretkklist[[i]][[which(pretkklist[[i]]=="FP"),2]]/(pretkklist[[i]][[which(pretkklist[[i]]=="FP"),2]]+0))
  } else {
    test.specs.t1[[i]][1,2] <- NA
  }
}

roc.pretkk <- rbindlist(test.specs.t1)
roc.pretkk <- roc.pretkk %>% mutate(time = "pre-T")%>%mutate_if(is.character, as.factor)

# 3wk
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
kk3wklist <- list()

# calculate how often each one occurs and get proportions 
for(i in 4:ncol(mergekk3wk)){
  kk3wklist[[i-3]] <-  mergekk3wk %>% group_by(mergekk3wk[,i]) %>% tally() %>%
    mutate(proportion = n/sum(n))
}

test.specs.kk <- list()
for(i in 1:length(kk3wklist)){
  test.specs.kk[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
  colnames(test.specs.kk[[i]]) <- c("TPR", "FPR")
  if("TP" %in% kk3wklist[[i]]$ `mergekk3wk[, i]` == TRUE && "FN" %in% kk3wklist[[i]]$ `mergekk3wk[, i]` == TRUE){
    test.specs.kk[[i]][1,1] <- kk3wklist[[i]][[which(kk3wklist[[i]]=="TP"),2]]/(kk3wklist[[i]][[which(kk3wklist[[i]]=="TP"),2]]+kk3wklist[[i]][[which(kk3wklist[[i]]=="FN"),2]])
  } else if("TP" %in% kk3wklist[[i]]$ `mergekk3wk[, i]`  == TRUE && "FN" %in% kk3wklist[[i]]$ `mergekk3wk[, i]` ==FALSE){
    test.specs.kk[[i]][1,1] <- kk3wklist[[i]][[which(kk3wklist[[i]]=="TP"),2]]/(kk3wklist[[i]][[which(kk3wklist[[i]]=="TP"),2]]+0)
  } else {
    test.specs.kk[[i]][1,1] <- NA
  }
  
  if("FP" %in% kk3wklist[[i]]$ `mergekk3wk[, i]` == TRUE && "TN" %in% kk3wklist[[i]]$ `mergekk3wk[, i]` ==TRUE){
    test.specs.kk[[i]][1,2] <- (kk3wklist[[i]][[which(kk3wklist[[i]]=="FP"),2]]/(kk3wklist[[i]][[which(kk3wklist[[i]]=="TN"),2]]+kk3wklist[[i]][[which(kk3wklist[[i]]=="FP"),2]]))
  } else if("FP" %in% kk3wklist[[i]]$ `mergekk3wk[, i]`  == TRUE && "TN" %in% kk3wklist[[i]]$ `mergekk3wk[, i]` ==FALSE){
    test.specs.kk[[i]][1,2] <- (kk3wklist[[i]][[which(kk3wklist[[i]]=="FP"),2]]/(kk3wklist[[i]][[which(kk3wklist[[i]]=="FP"),2]]+0))
  } else {
    test.specs.kk[[i]][1,2] <- NA
  }
}

roc.kk3wk <- rbindlist(test.specs.kk)
roc.kk3wk <- roc.kk3wk %>% mutate(time = "3 weeks")%>%mutate_if(is.character, as.factor)

# 9wk

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
kk9wklist <- list()

# calculate how often each one occurs and get proportions 
for(i in 4:ncol(mergekk9wk)){
  kk9wklist[[i-3]] <-  mergekk9wk %>% group_by(mergekk9wk[,i]) %>% tally() %>%
    mutate(proportion = n/sum(n))
}

test.specs.kk <- list()
for(i in 1:length(kk9wklist)){
  test.specs.kk[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
  colnames(test.specs.kk[[i]]) <- c("TPR", "FPR")
  if("TP" %in% kk9wklist[[i]]$ `mergekk9wk[, i]` == TRUE && "FN" %in% kk9wklist[[i]]$ `mergekk9wk[, i]` == TRUE){
    test.specs.kk[[i]][1,1] <- kk9wklist[[i]][[which(kk9wklist[[i]]=="TP"),2]]/(kk9wklist[[i]][[which(kk9wklist[[i]]=="TP"),2]]+kk9wklist[[i]][[which(kk9wklist[[i]]=="FN"),2]])
  } else if("TP" %in% kk9wklist[[i]]$ `mergekk9wk[, i]`  == TRUE && "FN" %in% kk9wklist[[i]]$ `mergekk9wk[, i]` ==FALSE){
    test.specs.kk[[i]][1,1] <- kk9wklist[[i]][[which(kk9wklist[[i]]=="TP"),2]]/(kk9wklist[[i]][[which(kk9wklist[[i]]=="TP"),2]]+0)
  } else {
    test.specs.kk[[i]][1,1] <- NA
  }
  
  if("FP" %in% kk9wklist[[i]]$ `mergekk9wk[, i]` == TRUE && "TN" %in% kk9wklist[[i]]$ `mergekk9wk[, i]` ==TRUE){
    test.specs.kk[[i]][1,2] <- (kk9wklist[[i]][[which(kk9wklist[[i]]=="FP"),2]]/(kk9wklist[[i]][[which(kk9wklist[[i]]=="TN"),2]]+kk9wklist[[i]][[which(kk9wklist[[i]]=="FP"),2]]))
  } else if("FP" %in% kk9wklist[[i]]$ `mergekk9wk[, i]`  == TRUE && "TN" %in% kk9wklist[[i]]$ `mergekk9wk[, i]` ==FALSE){
    test.specs.kk[[i]][1,2] <- (kk9wklist[[i]][[which(kk9wklist[[i]]=="FP"),2]]/(kk9wklist[[i]][[which(kk9wklist[[i]]=="FP"),2]]+0))
  } else {
    test.specs.kk[[i]][1,2] <- NA
  }
}

roc.kk9wk <- rbindlist(test.specs.kk)
roc.kk9wk <- roc.kk9wk %>% mutate(time = "9 weeks")%>%mutate_if(is.character, as.factor)

# 6m

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
kk6mlist <- list()

# calculate how often each one occurs and get proportions 
for(i in 4:ncol(mergekk6m)){
  kk6mlist[[i-3]] <-  mergekk6m %>% group_by(mergekk6m[,i]) %>% tally() %>%
    mutate(proportion = n/sum(n))
}

test.specs.kk <- list()
for(i in 1:length(kk6mlist)){
  test.specs.kk[[i]] <- as.data.frame(matrix(nrow=1, ncol=2))
  colnames(test.specs.kk[[i]]) <- c("TPR", "FPR")
  if("TP" %in% kk6mlist[[i]]$ `mergekk6m[, i]` == TRUE && "FN" %in% kk6mlist[[i]]$ `mergekk6m[, i]` == TRUE){
    test.specs.kk[[i]][1,1] <- kk6mlist[[i]][[which(kk6mlist[[i]]=="TP"),2]]/(kk6mlist[[i]][[which(kk6mlist[[i]]=="TP"),2]]+kk6mlist[[i]][[which(kk6mlist[[i]]=="FN"),2]])
  } else if("TP" %in% kk6mlist[[i]]$ `mergekk6m[, i]`  == TRUE && "FN" %in% kk6mlist[[i]]$ `mergekk6m[, i]` ==FALSE){
    test.specs.kk[[i]][1,1] <- kk6mlist[[i]][[which(kk6mlist[[i]]=="TP"),2]]/(kk6mlist[[i]][[which(kk6mlist[[i]]=="TP"),2]]+0)
  } else {
    test.specs.kk[[i]][1,1] <- NA
  }
  
  if("FP" %in% kk6mlist[[i]]$ `mergekk6m[, i]` == TRUE && "TN" %in% kk6mlist[[i]]$ `mergekk6m[, i]` ==TRUE){
    test.specs.kk[[i]][1,2] <- (kk6mlist[[i]][[which(kk6mlist[[i]]=="FP"),2]]/(kk6mlist[[i]][[which(kk6mlist[[i]]=="TN"),2]]+kk6mlist[[i]][[which(kk6mlist[[i]]=="FP"),2]]))
  } else if("FP" %in% kk6mlist[[i]]$ `mergekk6m[, i]`  == TRUE && "TN" %in% kk6mlist[[i]]$ `mergekk6m[, i]` ==FALSE){
    test.specs.kk[[i]][1,2] <- (kk6mlist[[i]][[which(kk6mlist[[i]]=="FP"),2]]/(kk6mlist[[i]][[which(kk6mlist[[i]]=="FP"),2]]+0))
  } else {
    test.specs.kk[[i]][1,2] <- NA
  }
}

roc.kk6m <- rbindlist(test.specs.kk)
roc.kk6m <- roc.kk6m %>% mutate(time = "6 months")%>%mutate_if(is.character, as.factor)

kk.spec <- bind_rows(roc.pretkk, roc.kk3wk, roc.kk9wk, roc.kk6m)%>% mutate_if(is.character, as.factor)
kk.spec$time <- factor(kk.spec$time, levels=c("pre-T", "3 weeks", "9 weeks", "6 months"))

kk.spec.qtiles <- kk.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs))

kk.spec.plot <- ggplot()+
  geom_jitter(kk.spec, mapping=aes(x=FPR, y=TPR, colour=time))+
  ylab("True Positive Rate") + xlab("False Positive Rate")+
  ggsave("kkspecplot.pdf")

#### CCA Model T positive ####
cca.status.sample <- status.kkcca[sample.int(nrow(status.kkcca), 500, replace=FALSE, prob=NULL),]

condition.data.frame <- mergekkcca %>% select(dateN, CID, ccanegpos)
ccanegpos.spec <- diag_spec2(cca.status.sample, condition.data.frame) 
specsensccaN <- ccanegpos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="negative")

condition.data.frame <- mergekkcca %>% select(dateN, CID, ccaTpos)
ccaTpos.spec <- diag_spec2(cca.status.sample, condition.data.frame) 
specsensccaT <- ccaTpos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="trace")

condition.data.frame <- mergekkcca %>% select(dateN, CID, cca1pos)
cca1pos.spec <- diag_spec2(cca.status.sample, condition.data.frame) 
specsenscca1 <- cca1pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs))%>% 
  mutate(threshold="+")

condition.data.frame <- mergekkcca %>% select(dateN, CID, cca2pos)
cca2pos.spec <- diag_spec2(cca.status.sample, condition.data.frame) 
specsenscca2 <- cca2pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs))%>% 
  mutate(threshold="++")

condition.data.frame <- mergekkcca %>% select(dateN, CID, cca3pos)
cca3pos.spec <- diag_spec2(cca.status.sample, condition.data.frame) 
specsenscca3 <- cca3pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs))%>% 
  mutate(threshold="+++")

ccasensspecquant <- bind_rows(specsensccaN, specsensccaT,specsenscca1, specsenscca2, specsenscca3 )%>%
  mutate_if(is.character, as.factor)%>%
  mutate(diagnostic="Plus(+)", threshold=factor(threshold, levels = c("negative", "trace", "+", "++", "+++")), 
         time=factor(time, levels=c("pre-T", "3 weeks", "9 weeks", "6 months")))


cca.spec.plot <- ggplot()+
  geom_ribbon(data=ccasensspecquant, aes(x=`FPR_50%`, ymin=`TPR_2.5%`, ymax=`TPR_97.5%`), fill="grey")+
  geom_jitter(ccasensspecquant, mapping=aes(x=`FPR_50%` , y=`TPR_50%`, colour=threshold))+
  facet_wrap(time~.)+
  ylab("True Positive Rate") + xlab("False Positive Rate")+
  theme_bw()+
  ggsave("roc.cca.pdf")

#### gscore ####

gscore.status.sample <- status.kkgscore[sample.int(nrow(status.kkgscore), 500, replace=FALSE, prob=NULL),]

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore1pos)
gscore1pos.spec <- diag_spec2(gscore.status.sample, condition.data.frame)
specsensgscore1 <- gscore1pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score1")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore2pos)
gscore2pos.spec <- diag_spec2(gscore.status.sample, condition.data.frame)
specsensgscore2 <- gscore2pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score2")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore3pos)
gscore3pos.spec <- diag_spec2(gscore.status.sample, condition.data.frame) 
specsensgscore3 <- gscore3pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score3")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore4pos)
gscore4pos.spec <- diag_spec2(gscore.status.sample, condition.data.frame)
specsensgscore4 <- gscore4pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score4")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore5pos)
gscore5pos.spec <- diag_spec2(gscore.status.sample, condition.data.frame)
specsensgscore5 <- gscore5pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score5")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore6pos)
gscore6pos.spec <- diag_spec2(gscore.status.sample, condition.data.frame) 
specsensgscore6 <- gscore6pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score6")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore7pos)
gscore7pos.spec <- diag_spec2(gscore.status.sample, condition.data.frame)
specsensgscore7 <- gscore7pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score7")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore8pos)
gscore8pos.spec <- diag_spec2(gscore.status.sample, condition.data.frame) 
specsensgscore8 <- gscore8pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score8")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore9pos)
gscore9pos.spec <- diag_spec2(gscore.status.sample, condition.data.frame)
specsensgscore9 <- gscore9pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score9")

condition.data.frame <- mergekkcca %>% select(dateN, CID, gscore10pos)
gscore10pos.spec <- diag_spec2(gscore.status.sample, condition.data.frame) 
specsensgscore10 <- gscore10pos.spec %>% group_by(time) %>%
  summarise_at(vars(TPR, FPR),funs(!!!p_funs)) %>% 
  mutate(threshold="G-Score10")

gscoresensspecquant <- bind_rows(specsensgscore1, specsensgscore2, specsensgscore3, specsensgscore4,
                                 specsensgscore5, specsensgscore6, specsensgscore7, specsensgscore8, 
                                 specsensgscore9, specsensgscore10)%>% 
                        mutate(time=factor(time, levels=c("pre-T", "3 weeks", "9 weeks", "6 months")))%>%
                        mutate(threshold=factor(threshold, levels=c("G-Score1", "G-Score2","G-Score3","G-Score4","G-Score5","G-Score6","G-Score7","G-Score8","G-Score9","G-Score10")))%>%
                        mutate(diagnostic="G-Score")

gscore.roc <- ggplot()+
  geom_ribbon(data=gscoresensspecquant, aes(x=`FPR_50%`, ymin=`TPR_2.5%`, ymax=`TPR_97.5%`),fill="grey", alpha=.6)+
  geom_jitter(gscoresensspecquant, mapping=aes(x=`FPR_50%` , y=`TPR_50%`,  colour=threshold))+
  facet_wrap(time~.)+
  ylab("Sensitivity") + xlab("False Positive Rate")+
  theme_bw()+
  ggsave("rocgscore.pdf")

cca.roc <- ggplot()+
  geom_ribbon(data=ccasensspecquant, aes(x=`FPR_50%`, ymin=`TPR_2.5%`, ymax=`TPR_97.5%`),fill="grey", alpha=.6)+
  geom_jitter(ccasensspecquant, mapping=aes(x=`FPR_50%` , y=`TPR_50%`,  colour=threshold))+
  facet_wrap(time~.)+
  ylab("Sensitivity") + xlab("False Positive Rate")+
  theme_bw()+
  ggsave("roccca.pdf")









