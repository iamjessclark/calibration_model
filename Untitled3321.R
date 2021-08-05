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
#source("posteriorfitting.R")

# get KK and CCA data
kkdata <- read.csv("KK.4tp.csv")
ccadata <- read.csv("cca.clean.csv")

kk <- loadKKdata("KK.4tp.csv")
KKIDs <- getKKChildIDs("KK.4tp.csv")

cca <- loadCCAdata("cca.clean.csv")
CCAIDs <- getCCAChildIDs("cca.clean.csv")
CID <- as.character(CCAIDs[(CCAIDs %in% KKIDs)])

kk <- kk[match(CID, KKIDs),,]
cca <- cca[match(CID, CCAIDs),]

# get G-Score data 
ccagscore <- loadgscoredata("cca.clean.csv")
gscoreCCAIDS <- getCCAChildIDsgscore("cca.clean.csv") 
gscoreCID <- as.character(gscoreCCAIDS[(gscoreCCAIDS %in% KKIDs)])
ccagscore <- ccagscore[match(gscoreCID,gscoreCCAIDS),]

#### Run models ####
#### KK and CCA Model ####

kkcca.model <- CalModKKCCA(210,4,6,kk,cca)

# monitored variables
prevKKcca1 <- density(c(kkcca.model$mcmc[[1]][,"prev[1]"], kkcca.model$mcmc[[2]][,"prev[1]"]))
prevKKcca2 <- density(c(kkcca.model$mcmc[[1]][,"prev[2]"], kkcca.model$mcmc[[2]][,"prev[2]"]))
prevKKcca3 <- density(c(kkcca.model$mcmc[[1]][,"prev[3]"], kkcca.model$mcmc[[2]][,"prev[3]"]))
prevKKcca4 <- density(c(kkcca.model$mcmc[[1]][,"prev[4]"], kkcca.model$mcmc[[2]][,"prev[4]"]))

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

#### KK and Gscore Model ####
kkGscore.model <- CalModKKGScore(210,4,6,kk,ccagscore)

prevKKgs1 <- density(c(kkGscore.model$mcmc[[1]][,"prev[1]"], kkGscore.model$mcmc[[2]][,"prev[1]"]))
prevKKgs2 <- density(c(kkGscore.model$mcmc[[1]][,"prev[2]"], kkGscore.model$mcmc[[2]][,"prev[2]"]))
prevKKgs3 <- density(c(kkGscore.model$mcmc[[1]][,"prev[3]"], kkGscore.model$mcmc[[2]][,"prev[3]"]))
prevKKgs4 <- density(c(kkGscore.model$mcmc[[1]][,"prev[4]"], kkGscore.model$mcmc[[2]][,"prev[4]"]))

rtnbKKgscore <- density(c(kkGscore.model$mcmc[[1]][,"rtnb"], kkGscore.model$mcmc[[2]][,"rtnb"]))

sht1KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkksh[1]"], kkGscore.model$mcmc[[2]][,"tkksh[1]"]))
sht2KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkksh[2]"], kkGscore.model$mcmc[[2]][,"tkksh[2]"]))
sht3KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkksh[3]"], kkGscore.model$mcmc[[2]][,"tkksh[3]"]))
sht4KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkksh[4]"], kkGscore.model$mcmc[[2]][,"tkksh[4]"]))

rtt1KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkkrt[1]"], kkGscore.model$mcmc[[2]][,"tkkrt[1]"]))
rtt2KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkkrt[2]"], kkGscore.model$mcmc[[2]][,"tkkrt[2]"]))
rtt3KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkkrt[3]"], kkGscore.model$mcmc[[2]][,"tkkrt[3]"]))
rtt4KKgscore <- density(c(kkGscore.model$mcmc[[1]][,"tkkrt[4]"], kkGscore.model$mcmc[[2]][,"tkkrt[4]"]))

klogKKgscore <- density(c(kkGscore.model$mcmc[[1]][,"k"], kkGscore.model$mcmc[[2]][,"k"]))
interceptKKgscore <- density(c(kkGscore.model$mcmc[[1]][,"intercept"], kkGscore.model$mcmc[[2]][,"intercept"]))

outputkkGscore <- as.data.frame(as.matrix(as.mcmc(kkGscore.model)))
status.kkgs <- outputkkGscore[,16:ncol(outputkkGscore)]

#### Figure 1 Raw Data ####

# kk data
kkpre <- getKKtime(kk,TimeStep=1)
kkthreewk <- getKKtime(kk, TimeStep = 2)
kkninewk <- getKKtime(kk, TimeStep = 3)
kksixmth <- getKKtime(kk, TimeStep = 4)

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
dplyr::select(CID, Gscore, CCA, kkstatus, meancount)%>%
  mutate(time="pre-treatment",
         ccanegpos = 1,
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
  dplyr::select(CID, Gscore, CCA, kkstatus, meancount)%>%
  mutate(time="three-weeks",
         ccanegpos = 1,
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
  dplyr::select(CID, Gscore, CCA, kkstatus, meancount)%>%
  mutate(time="nine-weeks",
         ccanegpos = 1,
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
  dplyr::select(CID, Gscore, CCA, kkstatus, meancount)%>%
  mutate(time="six-months",
         ccanegpos = 1,
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

cca_kk_gs <- bind_rows(predata, threewkdata, ninewkdata, sixmonthdata)%>%
  mutate(time = factor(time, levels=c("pre-treatment", "three-weeks", "nine-weeks","six-months")),
         epg=meancount*24,
         epg=ifelse(is.nan(epg), NA, epg),
         intensity=case_when(
           epg==0~"zero_count",
           epg>0 & epg <=99~"low",
           epg>=100 & epg <=399~"moderate",
           epg>=400 ~ "high", 
           TRUE~NA_character_),
         intensity = factor(intensity, levels=c("zero_count","low","moderate", "high")))

cca_kk_gs %>%
  filter(!is.na(CCA))%>%
  filter(!is.na(Gscore))%>%
  mutate(Gscore=as.factor(Gscore), CCA=as.factor(CCA),
    Gscore=fct_recode(Gscore, G1="0", G2="1", G3="2", G4="3", G5="4", G6="5", G7="6", G8="7", G9="8", G10="9"),
    CCA=fct_recode(CCA, Negative="0", Trace="1", "+"="2", "++"="3", "+++"="4"))%>%
  filter(!is.na(epg))%>%
  ggplot(aes(x=Gscore, y=CCA))+
  geom_jitter(aes(colour=intensity), size=2)+
  facet_grid(time~.)+
  labs(colour="Intensity Categories", y="CCA+", x="G-Score")+
  scale_colour_viridis(discrete=TRUE) +
  theme_bw()
  ggsave("fig1corr.pdf")

#### Calculating precision for normal draw of POC-CCAs ####

# cca and kk
ccaprec <- cca_kk_gs %>% 
  mutate(meancount, as.factor(meancount))%>%
  filter(!is.na(CCA))%>%
  filter(!is.nan(meancount))%>%
  group_by(meancount, CCA)%>%
  summarise(combos=n())%>%
  count(meancount)

kkccavar <- var(ccaprec$n)
sdevkkcca <- sqrt(kkccavar)# this is how many standard devations from the mean 
dev2 <- sdevkkcca/2 # divide by two to get how far on each side 
prec2 <- 1/(dev2^2)

# sens analysis

SAprec <- 1/(sdevkkcca^2)

# kk and G-Score
gsprec <- cca_kk_gs %>% 
  mutate(meancount, as.factor(meancount))%>%
  filter(!is.na(Gscore))%>%
  filter(!is.nan(meancount))%>%
  group_by(meancount, Gscore)%>%
  summarise(combos=n())%>%
  count(meancount)

count.var <- var(gsprec$n)
sdev <- sqrt(count.var)
dev <- sdev/2
prec <- 1/(dev^2)

# sens analysis
SAprecgs <- 1/(sdev^2)

#### Figure 2 Log Curve ####

# log curve

# cca sampling K and intercept 
kintkkcca <- as.data.frame(cbind(c(kkcca.model$mcmc[[1]][,"k"], kkcca.model$mcmc[[2]][,"k"]),c(kkcca.model$mcmc[[1]][,"intercept"], kkcca.model$mcmc[[2]][,"intercept"])))
kkcca.sample <- kintkkcca[sample.int(nrow(kintkkcca), 500, replace=FALSE, prob=NULL),]

# gscore sampling K and intercept 
kintkkgscore <- as.data.frame(cbind(c(kkGscore.model$mcmc[[1]][,"k"], kkGscore.model$mcmc[[2]][,"k"]),c(kkGscore.model$mcmc[[1]][,"intercept"], kkGscore.model$mcmc[[2]][,"intercept"])))
kkgscore.sample <- kintkkgscore[sample.int(nrow(kintkkgscore), 500, replace=FALSE, prob=NULL),]

logcurvelistkkcca <- list()
logcurvelistkkgscore <- list()
x=seq(from=1, to=100, by=1)

# cca
for(i in 1:100){
  logcurvelistkkcca[[i]] <- as.data.frame(logcurvecca(x,kintkkcca[i,1],kintkkcca[i,2]))
  logcurvelistkkcca[[i]] <- mutate(logcurvelistkkcca[[i]], id = rownames(logcurvelistkkcca[[i]]))
  colnames(logcurvelistkkcca[[i]]) <- c("logcurve", "epg")
}

logcca <- do.call(rbind, logcurvelistkkcca)

logcca$epg <- as.factor(logcca$epg)
logcca$sdmax <- logcca$logcurve+1.758821
logcca$sdmin <- logcca$logcurve-1.758821
ccaquants <- logcca %>% group_by(epg)%>%
  summarise(quantiles = quantile(logcurve, c(0.025, 0.5, 0.975)), 
            sdmax = quantile(sdmax, c(0.025, 0.5, 0.975)), 
            sdmin = quantile(sdmin, c(0.025, 0.5, 0.975)), 
            q = c(0.025, 0.5, 0.75)) %>%
  pivot_wider(names_from = q, values_from=c(quantiles, sdmax, sdmin))
ccaquants$sdmax_0.5[which(ccaquants$sdmax_0.5>4)] <- 4
ccaquants$epg <- as.numeric(as.character(ccaquants$epg))

# gscore
for(i in 1:100){
  logcurvelistkkgscore[[i]] <- as.data.frame(logcurvegscore(x,kintkkgscore[i,1],kintkkgscore[i,2]))
  logcurvelistkkgscore[[i]] <- mutate(logcurvelistkkgscore[[i]], id = rownames(logcurvelistkkgscore[[i]]))
  colnames(logcurvelistkkgscore[[i]]) <- c("logcurve", "epg")
}

loggscore <- do.call(rbind, logcurvelistkkgscore)

loggscore$epg <- as.factor(loggscore$epg)
loggscore$sdmax <- loggscore$logcurve+1.093606
loggscore$sdmin <- loggscore$logcurve-1.093606
gscorequants <- loggscore %>% group_by(epg)%>%
  summarise(quantiles = quantile(logcurve, c(0.025, 0.5, 0.975)), 
            sdmax = quantile(sdmax, c(0.025, 0.5, 0.975)), 
            sdmin = quantile(sdmin, c(0.025, 0.5, 0.975)), 
            q = c(0.025, 0.5, 0.75)) %>%
  pivot_wider(names_from = q, values_from=c(quantiles, sdmax, sdmin))
gscorequants$sdmax_0.5[which(gscorequants$sdmax_0.5>9)] <- 9
gscorequants$epg <- as.numeric(as.character(gscorequants$epg))

gscurve <- ggplot()+
  geom_line(data=gscorequants, aes(x=epg, y=quantiles_0.5), col="purple4")+
  geom_ribbon(data=gscorequants, aes(x=epg, ymin=sdmin_0.5, ymax=sdmax_0.5), fill="purple4", alpha=0.2)+
  theme_light()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_y_continuous(expand = expansion(mult = c(0, .02)), limits = c(0, 9), breaks=c(1:9))+
  scale_x_continuous(expand = expansion(mult = c(0, .02)), limits = c(0, 100))+
  ylab("")
  ggsave("logcurvegs.pdf")

ccacurve <- ggplot()+
  geom_line(data=ccaquants, aes(x=epg, y=quantiles_0.5), col="darkorange3")+
  geom_ribbon(data=ccaquants, aes(x=epg, ymin=sdmin_0.5, ymax=sdmax_0.5), fill="darkorange3", alpha=0.2)+
  theme_bw()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_y_continuous(expand = expansion(mult = c(0, .02)), limits = c(0, 4), breaks=c(1:9))+
  scale_x_continuous(expand = expansion(mult = c(0, .02)), limits = c(0, 100))+
  ylab("")
  ggsave("logcurvecca.pdf")

panel <- ggarrange(ccacurve,gscurve, labels=c("A", "B"), nrow=1, ncol=2)
ggsave("logcurvepanel.pdf", height=6, width=10)

 #### Figure 3 Prob Inf ####

ccainfprob <- time.steps(status.kkcca)
ccainfprob$diagnostic <- "CCA"
colnames(ccainfprob) <- c("time", "CID", "prob.inf", "diagnostic")
gsinfprob <- time.steps(status.kkgs)
gsinfprob$diagnostic <- "Gscore"
colnames(gsinfprob) <- c("time", "CID", "prob.inf", "diagnostic")

probsinf <- rbind(ccainfprob, gsinfprob)

zero <- cca_kk_gs %>% dplyr::select(time, Gscore, CCA, CID, intensity)%>%
  pivot_longer(!c(CID, time, intensity), names_to = "diagnostic", values_to = "score")%>%
  mutate(score=factor(score, levels=c(0,1,2,3,4,5,6,7,8,9)))%>%
  merge(probsinf, by=c("CID", "time", "diagnostic"))%>%
  dplyr::filter(intensity=="zero_count", !score==is.na(score))

col=c("#fdd0a2","#fd8d3c","#f16913","#d94801","#8c2d04")
ccaprobinfplot <- zero %>% filter(diagnostic=="CCA")%>%
  mutate(score = factor(score, levels=c("0", "1", "2", "3", "4"), 
                        labels=c("Negative", "Trace", "+", "++", "+++")))%>%
  ggplot(aes(x=prob.inf, y=score, colour=score))+
  geom_jitter(size=2, shape=1, stroke = 2)+
  labs(x="probability infected", y="POC-CCA+")+
  #scale_fill_manual(values=col)+
  scale_color_manual(values = col)+
  theme_light()+
  theme(legend.position = "left", text = element_text(size=16))
  ggsave("ccaxprobinf.pdf")

col=c("#dadaeb","#bcbddc","#bfd3e6","#9ebcda","#8c96c6","#8c6bb1","#88419d","#810f7c","#4d004b")
gsprobinf <- zero %>% filter(diagnostic=="Gscore")%>%
   mutate(score= as.numeric(as.character(score))+1,
     score = factor(score, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
                          labels=c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10")))%>%
  ggplot(aes(x=prob.inf, y=score, colour=score))+
  geom_jitter(size=2, shape=1, stroke = 2)+
  labs(x="probability infected", y="G-Score")+
  #scale_fill_manual(values=col)+
  scale_color_manual(values=col)+
  theme_light()+
  theme(text = element_text(size=16))
  ggsave("gscorexprobinf.pdf")

ggarrange(ccaprobinfplot, gsprobinf, labels=c("A", "B"), nrow=1, ncol=2)
ggsave("probinfpanel.pdf")


#### Figure 4 ROC Curves ####

# CCA #
cca.status.sample <- status.kkcca[sample.int(nrow(status.kkcca), 500, replace=FALSE, prob=NULL),]

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


# Pre-treatment 

predatacca <- predata %>% dplyr::select(CID, CCA) %>%
  right_join(t1cca, by="CID")

predatacca <- predatacca[-which(is.na(predatacca$CCA)==T),]
predictionspreT <- rep(list(predatacca$CCA), 500)
labelslistpreT <- list()
for(i in 3:ncol(predatacca)){
  labelslistpreT[[i-2]] <- predatacca[,i]
}
manypredccat1 <- prediction(predictionspreT, labelslistpreT)
auc.perfccat1 <- performance(manypredccat1, measure = "auc")
many.roc.ccat1 <-  performance(manypredccat1, measure = "tpr", x.measure = "fpr")

# three weeks 

threeweekscca <- threewkdata %>% dplyr::select(CID, CCA) %>%
  right_join(t2cca, by="CID")
threeweekscca <- threeweekscca[-which(is.na(threeweekscca$CCA)==T),]

predictionsthreewk <- rep(list(threeweekscca$CCA), 500)
labelslistthreewk <- list()
for(i in 3:ncol(threeweekscca)){
  labelslistthreewk[[i-2]] <- threeweekscca[,i]
}

manypredccat2 <- prediction(predictionsthreewk, labelslistthreewk)
auc.perfccat2 <- performance(manypredccat2, measure = "auc")
many.roc.ccat2 <-  performance(manypredccat2, measure = "tpr", x.measure = "fpr")

# nine weeks 

nineweekscca <- ninewkdata %>% dplyr::select(CID, CCA) %>%
  right_join(t3cca, by="CID")
nineweekscca <- nineweekscca[-which(is.na(nineweekscca$CCA)==T),]
predictionsninewk <- rep(list(nineweekscca$CCA), 500)
labelslistninewk <- list()
for(i in 3:ncol(nineweekscca)){
  labelslistninewk[[i-2]] <- nineweekscca[,i]
}

manypredccat3 <- prediction(predictionsninewk, labelslistninewk)
auc.perfccat3 <- performance(manypredccat3, measure = "auc")
many.roc.ccat3 <-  performance(manypredccat3, measure = "tpr", x.measure = "fpr")

# six months 

sixmonthscca <- sixmonthdata %>% dplyr::select(CID, CCA) %>%
  right_join(t4cca, by="CID")
sixmonthscca <- sixmonthscca[-which(is.na(sixmonthscca$CCA)==T),]
predictionssixmth <- rep(list(sixmonthscca$CCA), 500)
labelslistsixmonth <- list()
for(i in 3:ncol(sixmonthscca)){
  labelslistsixmonth[[i-2]] <- sixmonthscca[,i]
}

manypredccat4 <- prediction(predictionssixmth, labelslistsixmonth)
auc.perfccat4 <- performance(manypredccat4, measure = "auc")
many.roc.ccat4 <-  performance(manypredccat4, measure = "tpr", x.measure = "fpr")

ccaroct1 <- plotfunccca(many.roc.ccat1, auc.perfccat1)
ccaroct2 <- plotfunccca(many.roc.ccat2, auc.perfccat2)
ccaroct3 <- plotfunccca(many.roc.ccat3, auc.perfccat3)
ccaroct4 <- plotfunccca(many.roc.ccat4, auc.perfccat4)

# G-SCORE #

gs.status.sample <- status.kkgs[sample.int(nrow(status.kkgs), 500, replace=FALSE, prob=NULL),]

t1gs <- as.data.frame(gs.status.sample [,1:210])
t2gs <- as.data.frame(gs.status.sample[,211:420])
t3gs <- as.data.frame(gs.status.sample[,421:630])
t4gs <- as.data.frame(gs.status.sample[,631:840])

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


# Pre-treatment 

predatags <- predata %>% dplyr::select(CID, Gscore) %>%
  right_join(t1gs, by="CID")

predatags <- predatags[-which(is.na(predatags$Gscore)==T),]
predictionspreT <- rep(list(predatags$Gscore), 500)
labelslistpreT <- list()
for(i in 3:ncol(predatags)){
  labelslistpreT[[i-2]] <- predatags[,i]
}
manypredgst1 <- prediction(predictionspreT, labelslistpreT)
auc.perfgst1 <- performance(manypredgst1, measure = "auc")
many.roc.gst1 <-  performance(manypredgst1, measure = "tpr", x.measure = "fpr")

# three weeks 

threeweeksgs <- threewkdata %>% dplyr::select(CID, Gscore) %>%
  right_join(t2gs, by="CID")
threeweeksgs <- threeweeksgs[-which(is.na(threeweeksgs$Gscore)==T),]

predictionsthreewk <- rep(list(threeweeksgs$Gscore), 500)
labelslistthreewk <- list()
for(i in 3:ncol(threeweeksgs)){
  labelslistthreewk[[i-2]] <- threeweeksgs[,i]
}

manypredgst2 <- prediction(predictionsthreewk, labelslistthreewk)
auc.perfgst2 <- performance(manypredgst2, measure = "auc")
many.roc.gst2 <-  performance(manypredgst2, measure = "tpr", x.measure = "fpr")

# nine weeks 

nineweeksgs <- ninewkdata %>% dplyr::select(CID, Gscore) %>%
  right_join(t3gs, by="CID")
nineweeksgs <- nineweeksgs[-which(is.na(nineweeksgs$Gscore)==T),]
predictionsninewk <- rep(list(nineweeksgs$Gscore), 500)
labelslistninewk <- list()
for(i in 3:ncol(nineweeksgs)){
  labelslistninewk[[i-2]] <- nineweeksgs[,i]
}

manypredgst3 <- prediction(predictionsninewk, labelslistninewk)
auc.perfgst3 <- performance(manypredgst3, measure = "auc")
many.roc.gst3 <-  performance(manypredgst3, measure = "tpr", x.measure = "fpr")

# six months 

sixmonthsgs <- sixmonthdata %>% dplyr::select(CID, Gscore) %>%
  right_join(t4gs, by="CID")
sixmonthsgs <- sixmonthsgs[-which(is.na(sixmonthsgs$Gscore)==T),]
predictionssixmth <- rep(list(sixmonthsgs$Gscore), 500)
labelslistsixmonth <- list()
for(i in 3:ncol(sixmonthsgs)){
  labelslistsixmonth[[i-2]] <- sixmonthsgs[,i]
}

manypredgst4 <- prediction(predictionssixmth, labelslistsixmonth)
auc.perfgst4 <- performance(manypredgst4, measure = "auc")
many.roc.gst4 <-  performance(manypredgst4, measure = "tpr", x.measure = "fpr")

gsroct1 <- plotfuncgs(many.roc.gst1, auc.perfgst1)
gsroct2 <- plotfuncgs(many.roc.gst2, auc.perfgst2)
gsroct3 <- plotfuncgs(many.roc.gst3, auc.perfgst3)
gsroct4 <- plotfuncgst4(many.roc.gst4, auc.perfgst4)

ggarrange(ccaroct1,gsroct1, ccaroct2,gsroct2, ccaroct3,gsroct3, ccaroct4,gsroct4,
          ncol=2, 
          nrow=4, 
          labels = c("A","E", "B", "F", "C", "G","D", "H"))
  ggsave("roccurvespanel.pdf")

#### simulation ####
 
# functions for simulations 
source("simulations2.R")

prevalences <- runif(100, min=0.01, max=0.1)
meanval <- 1.66
LowListGS <- list()
nRep = 1
nFECrep = 2
nInd = 10000
nRuns=50
 
# G-Score 
for(i in 1:length(prevalences)){
 prev <- prevalences[i]
 LowListGS[[i]] <- t(sapply(1:nRuns,runRepeatsgs))
 print(i)
}

trueprevsgs <- extractprevs(LowListGS, prevalences, 1, 1, nInd, nRuns)
kkprevsgs <- extractkkprevs(LowListGS, prevalences, 2, 0, nInd, nRuns)
gscoreprevs <- extractprevs(LowListGS, prevalences,4,2, nInd, nRuns )

prevsgs <- as.data.frame(cbind(trueprevsgs,kkprevsgs, gscoreprevs))

prevsgs <- prevsgs[order(prevsgs$trueprevsgs),]
prevsgs$trueprevsgs <- as.factor(prevsgs$trueprevsgs)  

gsprevsgs <- prevsgs %>% 
  group_by(trueprevsgs) %>% 
  summarise(quantiles = quantile(gscoreprevs, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975))%>%
  pivot_wider(names_from = q, values_from = quantiles)%>%
  mutate(trueprevsgs=as.numeric(as.character(trueprevsgs)), diagnostic="G-Score")%>%
  rename(trueprev=trueprevsgs)

kkprevsgs <- prevsgs %>% 
  group_by(trueprevsgs) %>% 
  summarise(quantiles = quantile(kkprevsgs, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>%
  pivot_wider(names_from = q, values_from = quantiles)%>%
  mutate(trueprevsgs=as.numeric(as.character(trueprevsgs)), diagnostic="Kato-Katz" )%>%
  rename(trueprev=trueprevsgs)

simgsprevs <- bind_rows(gsprevsgs, kkprevsgs)

simgs <- ggplot()+
  geom_line(data=simgsprevs, aes(x=trueprevsgs*100, y=`0.5`*100, group=diagnostic, colour=diagnostic))+
  geom_ribbon(data=simgsprevs, aes(x=trueprevsgs*100, ymin=`0.025`*100, ymax=`0.975`*100, group=diagnostic, fill=diagnostic), alpha=0.2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red", alpha = 0.6)+
  theme_bw()+
  ylab("Estimated Prevalence %")+xlab("True Prevalence %")+
  scale_fill_manual(values=c("purple4","grey"))+
  scale_colour_manual(values=c("purple4","grey"))+
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 10))
  ggsave("gssimplot.pdf")

# CCA
LowListcca <- list()
for(i in 1:length(prevalences)){
  prev <- prevalences[i]
  LowListcca[[i]] <- t(sapply(1:nRuns,runRepeatscca))
  print(i)
}

trueprevscca <- extractprevs(LowListcca, prevalences, 1, 1, nInd, nRuns)
kkprevscca <- extractkkprevs(LowListcca, prevalences, 2, 0, nInd, nRuns)
ccaprevs <- extractprevs(LowListcca, prevalences,4, 2, nInd, nRuns )
ccaprevsTP <- extractprevs(LowListcca, prevalences,4, 1, nInd, nRuns )

prevscca <- as.data.frame(cbind(trueprevscca,kkprevscca, ccaprevs))

prevscca <- prevscca[order(prevscca$trueprevscca),]
prevscca$trueprevscca <- as.factor(prevscca$trueprevscca)  

ccaprevscca <- prevscca %>% 
  group_by(trueprevscca) %>% 
  summarise(quantiles = quantile(ccaprevs, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975))%>%
  pivot_wider(names_from = q, values_from = quantiles)%>%
  mutate(diagnostic="POC-CCA+" )%>%
  rename(trueprev=trueprevscca)

kkprevscca <- prevscca %>% 
  group_by(trueprevscca) %>% 
  summarise(quantiles = quantile(kkprevscca, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>%
  pivot_wider(names_from = q, values_from = quantiles)%>%
  mutate(trueprevscca=as.numeric(as.character(trueprevscca)), diagnostic="Kato-Katz" )%>%
  rename(trueprev=trueprevscca)

kkprevscca$trueprev <- as.factor(kkprevscca$trueprev )

prevsccaTP <- as.data.frame(cbind(trueprevscca,ccaprevsTP))
prevsccaTP <- prevsccaTP[order(prevsccaTP$trueprevscca),]
prevsccaTP$trueprevscca <- as.factor(prevsccaTP$trueprevscca) 
ccaprevsccaTP <- prevsccaTP %>% 
  group_by(trueprevscca) %>% 
  summarise(quantiles = quantile(ccaprevsTP, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975))%>%
  pivot_wider(names_from = q, values_from = quantiles)%>%
  mutate(trueprevscca=as.numeric(as.character(trueprevscca)), diagnostic="POC-CCA+ Trace Positive" )%>%
  rename(trueprev=trueprevscca)

simccaprevs <- bind_rows(ccaprevscca, kkprevscca)
simccaprevs$trueprev <- as.numeric(as.character(simccaprevs$trueprev))

simcca <- ggplot()+
  geom_line(data=simccaprevs, aes(x=trueprevscca*100, y=`0.5`*100, group=diagnostic, colour=diagnostic))+
  geom_ribbon(data=simccaprevs, aes(x=trueprevscca*100, ymin=`0.025`*100, ymax=`0.975`*100, group=diagnostic, fill=diagnostic), alpha=0.2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red", alpha = 0.6)+
  theme_bw()+
  ylab("Estimated Prevalence %")+xlab("True Prevalence %")+
  scale_fill_manual(values=c("grey", "darkorange3"))+
  scale_colour_manual(values=c("grey", "darkorange3"))+
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 10))
  ggsave("ccasimplot.pdf")

simpanel <- ggarrange(simcca, simgs, nrow=1, ncol=2, labels = c("A", "B"))
  ggsave("simpanel.pdf")

simgsprevs$trueprev <- as.factor(simgsprevs$trueprev )
simccaprevs$trueprev <- as.factor(simccaprevs$trueprev )

combineddata <- bind_rows(simccaprevs,simgsprevs )
combineddata$trueprev <- as.numeric(as.character(combineddata$trueprev))

simall <- ggplot()+
  geom_line(data=combineddata, aes(x=trueprev*100, y=`0.5`*100, group=diagnostic, colour=diagnostic))+
  geom_ribbon(data=combineddata, aes(x=trueprev*100, ymin=`0.025`*100, ymax=`0.975`*100, group=diagnostic, fill=diagnostic), alpha=0.2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red", alpha = 0.6)+
  theme_bw()+
  ylab("Estimated Prevalence %")+xlab("True Prevalence %")+
  scale_fill_manual(values=c("purple4","grey","darkorange3"))+
  scale_colour_manual(values=c("purple4","grey","darkorange3"))+
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 10))
  #ggsave("ccasimplot.pdf")

allsimdata <- bind_rows(simccaprevs, ccaprevsccaTP, simgsprevs)

ggplot()+
  geom_line(data=allsimdata, aes(x=trueprev*100, y=`0.5`*100, group=diagnostic, colour=diagnostic))+
  geom_ribbon(data=allsimdata, aes(x=trueprev*100, ymin=`0.025`*100, ymax=`0.975`*100, group=diagnostic, fill=diagnostic), alpha=0.2)+
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red", alpha = 0.6)+
  theme_bw()+
  ylab("Estimated Prevalence %")+xlab("True Prevalence %")+
  scale_fill_manual(values=c("purple4","grey","darkorange3", "cyan3"))+
  scale_colour_manual(values=c("purple4","grey","darkorange3", "cyan3"))+
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 10))
  ggsave("estorevssimplot.pdf")

#### simulated distribution of egg counts vs POC-CCA scores. ####
valslist <- getVals(LowListcca)
kkestsccamod <- melt(as.data.frame(valslist[[1]]))
colnames(kkestsccamod) <- c("variable", "kato-katz")
ccaestsccamod <- melt(as.data.frame(valslist[[2]]))
colnames(ccaestsccamod) <- c("variable", "cca")
ccamodsim <- cbind(kkestsccamod, ccaestsccamod)
ccamodsim <- ccamodsim[,-c(3)]
ccamodsim <- ccamodsim %>% 
  mutate(cca=as.factor(cca))
ccamodsim$`kato-katz` <- ccamodsim$`kato-katz`*24 #make into epg
ccamodsim$cca <- recode(ccamodsim$cca, `0`="Negative", `1`="Trace", `2`="+", `3`="++", `4`="+++")

ccamodsim <- ccamodsim %>% mutate(intensity=case_when(`kato-katz`==0 | `kato-katz`==0.0 ~ "zero count", 
                                                    `kato-katz`>0 & `kato-katz` <=99~"light infection", 
                                                    `kato-katz`>=100&`kato-katz`<=399 ~"moderate infection", 
                                                    TRUE~"heavy infection"), 
                                intensity=factor(intensity, levels=c("zero count", "light infection", "moderate infection", "heavy infection")))


propscca <- ccamodsim %>%
  group_by(cca, intensity)%>%
  summarise(n=n())%>%
  mutate(freq=100*n/sum(n))

cols=c("#d4b9da", "#c994c7", "#df65b0","#e7298a")
props1cca <- ggplot(data=propscca, aes(x=intensity, y=freq, fill=intensity)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position="bottom")+
  facet_grid(.~cca, switch="x")+
  scale_fill_manual(values=cols)+
  ylab("Proportion within score")+
  labs(fill = "WHO categories of \ninfection intensity")

props2cca <- ccamodsim %>%
  group_by(cca)%>%
  summarise(n=100*n()/nrow(.))

props2ccaplot <- ggplot(data=props2cca, aes(x=cca, y=n)) +
  geom_bar(stat="identity", position=position_dodge() , fill="#91003f")+
  theme_minimal()+
  theme(axis.title.x = element_blank())+
  ylab("Proportion of all scores")

simpanle2cca <- ggarrange(props2ccaplot,props1cca,  nrow=2, ncol=1, labels="auto", vjust=1, font.label = list(size=10))
ggsave("ccadistsimplotpanel.pdf")

# I don't know what this figure is #

colours=c( "#fdd0a2","#fd8d3c","#f16913","#d94801","#8c2d04",  "black")
lowsimcca <- ggplot(ccamodsim, aes(x = `kato-katz`, y = cca)) +
  geom_density_ridges(aes(fill=cca)) +
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16), 
        strip.text = element_text(size=12))+
  ylab("POC-CCA+") + xlab("eggs per gram of stool")+
  scale_fill_manual(values=colours)


# G-Score
gsvalslist <- getVals(LowListGS)
kkestsgsmod <- melt(as.data.frame(gsvalslist[[1]]))
colnames(kkestsgsmod) <- c("variable", "kato-katz")
gsestsgsmod <- as.data.frame(gsvalslist[[2]])
gsestsgsmod <- melt(gsestsgsmod)
colnames(gsestsgsmod) <- c("variable", "gs")
gsmodsim <- cbind(kkestsgsmod, gsestsgsmod)
gsmodsim <- gsmodsim[,-c(3)]
gsmodsim <- gsmodsim %>% 
  mutate(gs=as.factor(gs))
gsmodsim$`kato-katz` <- gsmodsim$`kato-katz`*24 #make into epg
gsmodsim$gs <- recode(gsmodsim$gs, `0`="G1", `1`="G2", `2`="G3", `3`="G4", `4`="G5", `5`="G6",`6`="G7", `7`="G8", `8`="G9", `9`="G10")
gsmodsim <- gsmodsim %>% mutate(intensity=case_when(`kato-katz`==0 | `kato-katz`==0.0 ~ "zero count", 
                             `kato-katz`>0 & `kato-katz` <=99~"light infection", 
                             `kato-katz`>=100&`kato-katz`<=399 ~"moderate infection", 
                             TRUE~"heavy infection"), 
                             intensity=factor(intensity, levels=c("zero count", "light infection", "moderate infection", "heavy infection")))
props <- gsmodsim %>%
  group_by(gs, intensity)%>%
  summarise(n=n())%>%
  mutate(freq=100*n/sum(n))


cols=c("#99d8c9","#66c2a4","#2ca25f","#006d2c")
props1 <- ggplot(data=props, aes(x=intensity, y=freq, fill=intensity)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position="bottom")+
  facet_grid(.~gs, switch="x")+
  scale_fill_manual(values=cols)+
  ylab("Proportion within score")+
  labs(fill = "WHO categories of \ninfection intensity")

props2 <- gsmodsim %>%
  group_by(gs)%>%
  summarise(n=100*n()/nrow(.))

props2plot <- ggplot(data=props2, aes(x=gs, y=n)) +
  geom_bar(stat="identity", position=position_dodge() , fill="#1c9099")+
  theme_minimal()+
  theme(axis.title.x = element_blank())+
  ylab("Proportion of all scores")


colours=c("#dadaeb","#bcbddc","#bfd3e6","#9ebcda","#74a9cf", "#8c96c6","#8c6bb1","#88419d","#810f7c", "black")
lowsimgs <- ggplot(gsmodsim, aes(x = `kato-katz`, y = gs)) +
  geom_density_ridges(aes(fill=gs))+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16), 
        strip.text = element_text(size=12))+
  ylab("G-Score+") + xlab("eggs per gram of stool")+
  scale_fill_manual(values=colours)



distlowpanel <- ggarrange(lowsimcca, lowsimgs, nrow=1, ncol=2, labels=c("A", "B"))
ggsave("distlowpanel.pdf")

simpanelall <- ggarrange(simpanel, distlowpanel, nrow=2, ncol=1)
ggsave("simallpanel.pdf")

gssimpanle2 <- ggarrange(props2plot,props1, nrow=2, ncol=1, labels="auto", vjust=1, font.label = list(size=10))
ggsave("distsimplotpanel.pdf")

#### Sensitivity Analysis #### 


#### Posteriors ####

# prevalences 

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="Prevalence CCA+", ylab="Scaled Density", cex.axis=1)
lines(prevKKcca1$x,prevKKcca1$y/max(prevKKcca1$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(prevKKcca1$x), prevKKcca1$x), c(rev(prevKKcca1$y/max(prevKKcca1$y)), rep(0,length(prevKKcca1$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(prevKKcca2$x,prevKKcca2$y/max(prevKKcca2$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(prevKKcca2$x), prevKKcca2$x), c(rev(prevKKcca2$y/max(prevKKcca2$y)), rep(0,length(prevKKcca2$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(prevKKcca3$x,prevKKcca3$y/max(prevKKcca3$y), lwd = 2.5, col="thistle1")
polygon(c(rev(prevKKcca3$x), prevKKcca3$x), c(rev(prevKKcca3$y/max(prevKKcca3$y)), rep(0,length(prevKKcca3$y))),
        col = adjustcolor("thistle1", alpha=.3))
lines(prevKKcca4$x,prevKKcca4$y/max(prevKKcca4$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(prevKKcca4$x), prevKKcca4$x), c(rev(prevKKcca4$y/max(prevKKcca4$y)), rep(0,length(prevKKcca4$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topleft", c("preT", "3 week", "9 week", "6 months"),
       col=c("#1b9e77","#d95f02","thistle1", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "prevsccaptp.pdf", height = 6, width = 6)
dev.off()
graphics.off()


dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="Prevalence G-Score", ylab="Scaled Density", cex.axis=1)
lines(prevKKgs1$x,prevKKgs1$y/max(prevKKgs1$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(prevKKgs1$x), prevKKgs1$x), c(rev(prevKKgs1$y/max(prevKKgs1$y)), rep(0,length(prevKKgs1$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(prevKKgs2$x,prevKKgs2$y/max(prevKKgs2$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(prevKKgs2$x), prevKKgs2$x), c(rev(prevKKgs2$y/max(prevKKgs2$y)), rep(0,length(prevKKgs2$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(prevKKgs3$x,prevKKgs3$y/max(prevKKgs3$y), lwd = 2.5, col="thistle1")
polygon(c(rev(prevKKgs3$x), prevKKgs3$x), c(rev(prevKKgs3$y/max(prevKKgs3$y)), rep(0,length(prevKKgs3$y))),
        col = adjustcolor("thistle1", alpha=.3))
lines(prevKKgs4$x,prevKKgs4$y/max(prevKKgs4$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(prevKKgs4$x), prevKKgs4$x), c(rev(prevKKgs4$y/max(prevKKgs4$y)), rep(0,length(prevKKgs4$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topleft", c("preT", "3 week", "9 week", "6 months"),
       col=c("#1b9e77","#d95f02","thistle1", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "prevsgsptp.pdf", height = 6, width = 6)
dev.off()
graphics.off()



#K parameters 

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="Growth Parameter k", ylab="Scaled Density", cex.axis=1)
lines(klogKKcca$x,klogKKcca$y/max(klogKKcca$y), lwd = 2.5, col="darkorange3")
polygon(c(rev(klogKKcca$x), klogKKcca$x), c(rev(klogKKcca$y/max(klogKKcca$y)), rep(0,length(klogKKcca$y))),
        col = adjustcolor("darkorange3", alpha=.3))
lines(klogKKgscore$x,klogKKgscore$y/max(klogKKgscore$y), lwd = 2.5, col="purple4")
polygon(c(rev(klogKKgscore$x), klogKKgscore$x), c(rev(klogKKgscore$y/max(klogKKgscore$y)), rep(0,length(klogKKgscore$y))),
        col = adjustcolor("purple4", alpha=.3))

legend("topleft", c("Kato-Katz & CCA+", "KK & G-Score"),
       col=c("darkorange3","purple4"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "kparsposts.pdf", height = 6, width = 6)
dev.off()
graphics.off()


# intercepts

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(-1.5,3), ylim = c(0,1), xlab="Intercept", ylab="Scaled Density", cex.axis=1)
lines(interceptKKcca$x,interceptKKcca$y/max(interceptKKcca$y), lwd = 2.5, col="darkorange3")
polygon(c(rev(interceptKKcca$x), interceptKKcca$x), c(rev(interceptKKcca$y/max(interceptKKcca$y)), rep(0,length(interceptKKcca$y))),
        col = adjustcolor("darkorange3", alpha=.3))
lines(interceptKKgscore$x,interceptKKgscore$y/max(interceptKKgscore$y), lwd = 2.5, col="purple4")
polygon(c(rev(interceptKKgscore$x), interceptKKgscore$x), c(rev(interceptKKgscore$y/max(interceptKKgscore$y)), rep(0,length(interceptKKgscore$y))),
        col = adjustcolor("purple4", alpha=.3))

legend("topleft", c("Kato-Katz & CCA+", "KK & G-Score"),
       col=c("darkorange3","purple4"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "intparsposts.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# shapes cca

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="shape parameters CCA+", ylab="Scaled Density", cex.axis=1)
lines(sht1KKcca$x,sht1KKcca$y/max(sht1KKcca$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(sht1KKcca$x), sht1KKcca$x), c(rev(sht1KKcca$y/max(sht1KKcca$y)), rep(0,length(sht1KKcca$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(sht2KKcca$x,sht2KKcca$y/max(sht2KKcca$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(sht2KKcca$x), sht2KKcca$x), c(rev(sht2KKcca$y/max(sht2KKcca$y)), rep(0,length(sht2KKcca$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(sht3KKcca$x,sht3KKcca$y/max(sht3KKcca$y), lwd = 2.5, col="thistle1")
polygon(c(rev(sht3KKcca$x), sht3KKcca$x), c(rev(sht3KKcca$y/max(sht3KKcca$y)), rep(0,length(sht3KKcca$y))),
        col = adjustcolor("thistle1", alpha=.3))
lines(sht4KKcca$x,sht4KKcca$y/max(sht4KKcca$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(sht4KKcca$x), sht4KKcca$x), c(rev(sht4KKcca$y/max(sht4KKcca$y)), rep(0,length(sht4KKcca$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topright", c("preT", "3 week", "9 week", "6 months"),
       col=c("#1b9e77","#d95f02","thistle1", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "shapeccaptp.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# shapes gscore

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="shape parameters G-Score", ylab="Scaled Density", cex.axis=1)
lines(sht1KKgscore$x,sht1KKgscore$y/max(sht1KKgscore$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(sht1KKgscore$x), sht1KKgscore$x), c(rev(sht1KKgscore$y/max(sht1KKgscore$y)), rep(0,length(sht1KKgscore$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(sht2KKgscore$x,sht2KKgscore$y/max(sht2KKgscore$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(sht2KKgscore$x), sht2KKgscore$x), c(rev(sht2KKgscore$y/max(sht2KKgscore$y)), rep(0,length(sht2KKgscore$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(sht3KKgscore$x,sht3KKgscore$y/max(sht3KKgscore$y), lwd = 2.5, col="thistle1")
polygon(c(rev(sht3KKgscore$x), sht3KKgscore$x), c(rev(sht3KKgscore$y/max(sht3KKgscore$y)), rep(0,length(sht3KKgscore$y))),
        col = adjustcolor("thistle1", alpha=.3))
lines(sht4KKgscore$x,sht4KKgscore$y/max(sht4KKgscore$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(sht4KKgscore$x), sht4KKgscore$x), c(rev(sht4KKgscore$y/max(sht4KKgscore$y)), rep(0,length(sht4KKgscore$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topright", c("preT", "3 week", "9 week", "6 months"),
       col=c("#1b9e77","#d95f02","thistle1", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "shapegsptp.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# rates cca

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="rate parameters CCA+", ylab="Scaled Density", cex.axis=1)
lines(rtt1KKcca$x,rtt1KKcca$y/max(rtt1KKcca$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(rtt1KKcca$x), rtt1KKcca$x), c(rev(rtt1KKcca$y/max(rtt1KKcca$y)), rep(0,length(rtt1KKcca$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(rtt2KKcca$x,rtt2KKcca$y/max(rtt2KKcca$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(rtt2KKcca$x), rtt2KKcca$x), c(rev(rtt2KKcca$y/max(rtt2KKcca$y)), rep(0,length(rtt2KKcca$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(rtt3KKcca$x,rtt3KKcca$y/max(rtt3KKcca$y), lwd = 2.5, col="thistle1")
polygon(c(rev(rtt3KKcca$x), rtt3KKcca$x), c(rev(rtt3KKcca$y/max(rtt3KKcca$y)), rep(0,length(rtt3KKcca$y))),
        col = adjustcolor("thistle1", alpha=.3))
lines(rtt4KKcca$x,rtt4KKcca$y/max(rtt4KKcca$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(rtt4KKcca$x), rtt4KKcca$x), c(rev(rtt4KKcca$y/max(rtt4KKcca$y)), rep(0,length(rtt4KKcca$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topright", c("preT", "3 week", "9 week", "6 months"),
       col=c("#1b9e77","#d95f02","thistle1", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "rateccaptp.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# rates gscore

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="rate parameters G-Score", ylab="Scaled Density", cex.axis=1)
lines(rtt1KKgscore$x,rtt1KKgscore$y/max(rtt1KKgscore$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(rtt1KKgscore$x), rtt1KKgscore$x), c(rev(rtt1KKgscore$y/max(rtt1KKgscore$y)), rep(0,length(rtt1KKgscore$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(rtt2KKgscore$x,rtt2KKgscore$y/max(rtt2KKgscore$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(rtt2KKgscore$x), rtt2KKgscore$x), c(rev(rtt2KKgscore$y/max(rtt2KKgscore$y)), rep(0,length(rtt2KKgscore$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(rtt3KKgscore$x,rtt3KKgscore$y/max(rtt3KKgscore$y), lwd = 2.5, col="thistle1")
polygon(c(rev(rtt3KKgscore$x), rtt3KKgscore$x), c(rev(rtt3KKgscore$y/max(rtt3KKgscore$y)), rep(0,length(rtt3KKgscore$y))),
        col = adjustcolor("thistle1", alpha=.3))
lines(rtt4KKgscore$x,rtt4KKgscore$y/max(rtt4KKgscore$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(rtt4KKgscore$x), rtt4KKgscore$x), c(rev(rtt4KKgscore$y/max(rtt4KKgscore$y)), rep(0,length(rtt4KKgscore$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topright", c("preT", "3 week", "9 week", "6 months"),
       col=c("#1b9e77","#d95f02","thistle1", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "rategsptp.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# negative binomial rate parameter 

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1.2), ylim = c(0,1), xlab="Negative Binomial Shape Parameter", ylab="Scaled Density", cex.axis=1)
lines(rtnbKKcca$x,rtnbKKcca$y/max(rtnbKKcca$y), lwd = 2.5, col="darkorange3")
polygon(c(rev(rtnbKKcca$x), rtnbKKcca$x), c(rev(rtnbKKcca$y/max(rtnbKKcca$y)), rep(0,length(rtnbKKcca$y))),
        col = adjustcolor("darkorange3", alpha=.3))
lines(rtnbKKgscore$x,rtnbKKgscore$y/max(rtnbKKgscore$y), lwd = 2.5, col="purple4")
polygon(c(rev(rtnbKKgscore$x), rtnbKKgscore$x), c(rev(rtnbKKgscore$y/max(rtnbKKgscore$y)), rep(0,length(rtnbKKgscore$y))),
        col = adjustcolor("purple4", alpha=.3))

legend("topleft", c("Kato-Katz & CCA+", "KK & G-Score"),
       col=c("darkorange3","purple4"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "rtnbposts.pdf", height = 6, width = 6)
dev.off()
graphics.off()

#save.image(file="diagnostic_model_environment.RData")
