SensAnCCA <- CalModKKCCA(210,4,6,kk,cca)
ccaoutput <- as.data.frame(as.matrix(as.mcmc(SensAnCCA)))

preccca1 <- density(c(SensAnCCA$mcmc[[1]][,"est_prec[1]"], SensAnCCA$mcmc[[2]][,"est_prec[1]"]))
preccca2 <- density(c(SensAnCCA$mcmc[[1]][,"est_prec[2]"], SensAnCCA$mcmc[[2]][,"est_prec[2]"]))
preccca3 <- density(c(SensAnCCA$mcmc[[1]][,"est_prec[3]"], SensAnCCA$mcmc[[2]][,"est_prec[3]"]))
preccca4 <- density(c(SensAnCCA$mcmc[[1]][,"est_prec[4]"], SensAnCCA$mcmc[[2]][,"est_prec[4]"]))

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,10), ylim = c(0,1), xlab="estimated precision parameter", ylab="Scaled Density", cex.axis=1)
lines(preccca1$x,preccca1$y/max(preccca1$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(preccca1$x), preccca1$x), c(rev(preccca1$y/max(preccca1$y)), rep(0,length(preccca1$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(preccca2$x,preccca2$y/max(preccca2$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(preccca2$x), preccca2$x), c(rev(preccca2$y/max(preccca2$y)), rep(0,length(preccca2$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(preccca3$x,preccca3$y/max(preccca3$y), lwd = 2.5, col="thistle1")
polygon(c(rev(preccca3$x), preccca3$x), c(rev(preccca3$y/max(preccca3$y)), rep(0,length(preccca3$y))),
        col = adjustcolor("thistle1", alpha=.3))
lines(preccca4$x,preccca4$y/max(preccca4$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(preccca4$x), preccca4$x), c(rev(preccca4$y/max(preccca4$y)), rep(0,length(preccca4$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topright", c("preT", "3 week", "9 week", "6 months"),
       col=c("#1b9e77","#d95f02","thistle1", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "sensanCCA.pdf", height = 6, width = 6)
dev.off()
graphics.off()

#### 2 standard deviations ####

CCA2SDs <- CalModKKCCA(210,4,6,kk,cca)
SD2ccaoutput <- as.data.frame(as.matrix(as.mcmc(CCA2SDs)))

prevKKccasa <- density(c(CCA2SDs$mcmc[[1]][,"prev[1]"], CCA2SDs$mcmc[[2]][,"prev[1]"]))
prevKKccasa <- density(c(CCA2SDs$mcmc[[1]][,"prev[2]"], CCA2SDs$mcmc[[2]][,"prev[2]"]))
prevKKccasa <- density(c(CCA2SDs$mcmc[[1]][,"prev[3]"], CCA2SDs$mcmc[[2]][,"prev[3]"]))
prevKKccasa <- density(c(CCA2SDs$mcmc[[1]][,"prev[4]"], CCA2SDs$mcmc[[2]][,"prev[4]"]))

rtnbKKccasa <- density(c(CCA2SDs$mcmc[[1]][,"rtnb"], CCA2SDs$mcmc[[2]][,"rtnb"]))
klogKKccasa <- density(c(CCA2SDs$mcmc[[1]][,"k"], CCA2SDs$mcmc[[2]][,"k"]))
interceptKKccasa <- density(c(CCA2SDs$mcmc[[1]][,"intercept"], CCA2SDs$mcmc[[2]][,"intercept"]))

sht1KKccasa <- density(c(CCA2SDs$mcmc[[1]][,"tkksh[1]"], CCA2SDs$mcmc[[2]][,"tkksh[1]"]))
sht2KKccasa <- density(c(CCA2SDs$mcmc[[1]][,"tkksh[2]"], CCA2SDs$mcmc[[2]][,"tkksh[2]"]))
sht3KKccasa <- density(c(CCA2SDs$mcmc[[1]][,"tkksh[3]"], CCA2SDs$mcmc[[2]][,"tkksh[3]"]))
sht4KKccasa <- density(c(CCA2SDs$mcmc[[1]][,"tkksh[4]"], CCA2SDs$mcmc[[2]][,"tkksh[4]"]))

rtt1KKccasa <- density(c(CCA2SDs$mcmc[[1]][,"tkkrt[1]"], CCA2SDs$mcmc[[2]][,"tkkrt[1]"]))
rtt2KKccasa <- density(c(CCA2SDs$mcmc[[1]][,"tkkrt[2]"], CCA2SDs$mcmc[[2]][,"tkkrt[2]"]))
rtt3KKccasa <- density(c(CCA2SDs$mcmc[[1]][,"tkkrt[3]"], CCA2SDs$mcmc[[2]][,"tkkrt[3]"]))
rtt4KKccasa <- density(c(CCA2SDs$mcmc[[1]][,"tkkrt[4]"], CCA2SDs$mcmc[[2]][,"tkkrt[4]"]))
status.kkccasa <- SD2ccaoutput[,16:ncol(SD2ccaoutput)]

#K parameters 

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="Prevalence", ylab="Scaled Density", cex.axis=1)
lines(klogKKccasa$x,klogKKccasa$y/max(klogKKccasa$y), lwd = 2.5, col="darkorange3")
polygon(c(rev(klogKKccasa$x), klogKKccasa$x), c(rev(klogKKccasa$y/max(klogKKccasa$y)), rep(0,length(klogKKccasa$y))),
        col = adjustcolor("darkorange3", alpha=.3))

legend("topleft", c("Kato-Katz & CCA+"),
       col=c("darkorange3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "SAkparsposts.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# intercept
dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(-2,1), ylim = c(0,1), xlab="Prevalence", ylab="Scaled Density", cex.axis=1)
lines(interceptKKccasa$x,interceptKKccasa$y/max(interceptKKccasa$y), lwd = 2.5, col="darkorange3")
polygon(c(rev(interceptKKccasa$x), interceptKKccasa$x), c(rev(interceptKKccasa$y/max(interceptKKccasa$y)), rep(0,length(interceptKKccasa$y))),
        col = adjustcolor("darkorange3", alpha=.3))


legend("topleft", c("Kato-Katz & CCA+"),
       col=c("darkorange3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "SAintparsposts.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# shapes 

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="shape parameters CCA+", ylab="Scaled Density", cex.axis=1)
lines(sht1KKccasa$x,sht1KKccasa$y/max(sht1KKccasa$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(sht1KKccasa$x), sht1KKccasa$x), c(rev(sht1KKccasa$y/max(sht1KKccasa$y)), rep(0,length(sht1KKccasa$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(sht2KKccasa$x,sht2KKccasa$y/max(sht2KKccasa$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(sht2KKccasa$x), sht2KKccasa$x), c(rev(sht2KKccasa$y/max(sht2KKccasa$y)), rep(0,length(sht2KKccasa$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(sht3KKccasa$x,sht3KKccasa$y/max(sht3KKccasa$y), lwd = 2.5, col="thistle1")
polygon(c(rev(sht3KKccasa$x), sht3KKccasa$x), c(rev(sht3KKccasa$y/max(sht3KKccasa$y)), rep(0,length(sht3KKccasa$y))),
        col = adjustcolor("thistle1", alpha=.3))
lines(sht4KKccasa$x,sht4KKccasa$y/max(sht4KKccasa$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(sht4KKccasa$x), sht4KKccasa$x), c(rev(sht4KKccasa$y/max(sht4KKccasa$y)), rep(0,length(sht4KKccasa$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topright", c("preT", "3 week", "9 week", "6 months"),
       col=c("#1b9e77","#d95f02","thistle1", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "SAshapeccaptp.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# rates

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1), ylim = c(0,1), xlab="rate parameters CCA+", ylab="Scaled Density", cex.axis=1)
lines(rtt1KKccasa$x,rtt1KKccasa$y/max(rtt1KKccasa$y), lwd = 2.5, col="#1b9e77")
polygon(c(rev(rtt1KKccasa$x), rtt1KKccasa$x), c(rev(rtt1KKccasa$y/max(rtt1KKccasa$y)), rep(0,length(rtt1KKccasa$y))),
        col = adjustcolor("#1b9e77", alpha=.3))
lines(rtt2KKccasa$x,rtt2KKccasa$y/max(rtt2KKccasa$y), lwd = 2.5, col="#d95f02")
polygon(c(rev(rtt2KKccasa$x), rtt2KKccasa$x), c(rev(rtt2KKccasa$y/max(rtt2KKccasa$y)), rep(0,length(rtt2KKccasa$y))),
        col = adjustcolor("#d95f02", alpha=.3))
lines(rtt3KKccasa$x,rtt3KKccasa$y/max(rtt3KKccasa$y), lwd = 2.5, col="thistle1")
polygon(c(rev(rtt3KKccasa$x), rtt3KKccasa$x), c(rev(rtt3KKccasa$y/max(rtt3KKccasa$y)), rep(0,length(rtt3KKccasa$y))),
        col = adjustcolor("thistle1", alpha=.3))
lines(rtt4KKccasa$x,rtt4KKccasa$y/max(rtt4KKccasa$y), lwd = 2.5, col="#7570b3")
polygon(c(rev(rtt4KKccasa$x), rtt4KKccasa$x), c(rev(rtt4KKccasa$y/max(rtt4KKccasa$y)), rep(0,length(rtt4KKccasa$y))),
        col = adjustcolor("#7570b3", alpha=.3))
legend("topright", c("preT", "3 week", "9 week", "6 months"),
       col=c("#1b9e77","#d95f02","thistle1", "#7570b3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "SArateccaptp.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# nb rate

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(0,1.2), ylim = c(0,1), xlab="Shape Negative Binomial Parameter", ylab="Scaled Density", cex.axis=1)
lines(rtnbKKccasa$x,rtnbKKccasa$y/max(rtnbKKccasa$y), lwd = 2.5, col="darkorange3")
polygon(c(rev(rtnbKKccasa$x), rtnbKKccasa$x), c(rev(rtnbKKccasa$y/max(rtnbKKccasa$y)), rep(0,length(rtnbKKccasa$y))),
        col = adjustcolor("darkorange3", alpha=.3))

legend("topleft", c("Kato-Katz & CCA+"),
       col=c("darkorange3"), 
       lty=c(1), lwd=2.5, bty='n',cex=0.75)

dev.copy(pdf, "SArtnbposts.pdf", height = 6, width = 6)
dev.off()
graphics.off()

# prob inf 

saccainfprob <- time.steps(status.kkccasa)
saccainfprob$diagnostic <- "CCA"
colnames(saccainfprob) <- c("time", "CID", "prob.inf", "diagnostic")

zero <- cca_kk_gs %>% dplyr::select(time, Gscore, CCA, CID, intensity)%>%
  pivot_longer(!c(CID, time, intensity), names_to = "diagnostic", values_to = "score")%>%
  mutate(score=factor(score, levels=c(0,1,2,3,4,5,6,7,8,9)))%>%
  merge(saccainfprob, by=c("CID", "time", "diagnostic"))%>%
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
ggsave("SAccaxprobinf.pdf")
