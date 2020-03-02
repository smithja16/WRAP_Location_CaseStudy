#Code for Location WRAP Workshop
#Compare models using simulated data

#JS: *** this version explores issues related to poor fit of warmer part of temperature response curve

#----Directories----
#Set your working directory!
# setwd('~/PROJECTS/WRAP Location/')

dir.create('Sim1') #data and model outputs will be saved locally here
Sim1 <- ('Sim1/')

#----Load Library & Function----
library(dismo)
library(gbm)
library(mgcv)
library(ggplot2)
library(viridis)
# source('WRAP_Location_CaseStudy/SimulatedWorld_Function.R') #load simulation function

#-----Simulate data----
dat <- SimulateWorld(temp_diff = 2, 
                     temp_spatial = "matern", 
                     PA_shape = "linear", 
                     abund_enviro = "poisson")  #takes a few minutes
#SB's original: temp_diff = 4, temp_spatial = "simple", PA_shape = "logistic", abund_enviro = "lnorm_low"
#EW's update: temp_diff = 4, temp_spatial = "matern", PA_shape = "logistic", abund_enviro = "lnorm_high"
#JS's sim: temp_diff = 2, temp_spatial = "matern", PA_shape = "linear", abund_enviro = "poisson"

colnames(dat)[1:2] <- c("Lon","Lat")
# saveRDS(dat, paste0(Sim1,'Sim1.rds')) #save data 
# # dat <- readRDS(paste0(Sim1,'Sim1.rds')) #load in data if needed

#Create dataframe with historical/forecast data
dat_hist <- dat[dat$year<=2020,]
dat_fcast <- dat[dat$year>2020,]

# #Make some quick plots to explore the data
# #All Years
# par(mfrow=c(2,2))
# plot(aggregate(suitability~year,dat,FUN="mean"),type="l", lwd=2, ylab="Suitability",col="dark grey", ylim=c(0,1))
# lines(aggregate(suitability~year,dat[dat$year<=2020,],FUN="mean"),col="blue")
# lines(aggregate(suitability~year,dat[dat$year<=2020,],FUN="min"),col="red")  #JS****
# plot(aggregate(pres~year,dat,FUN="mean"),type="l", lwd=2,ylab="Presence",col="dark grey", ylim=c(0,1))
# lines(aggregate(pres~year,dat[dat$year<=2020,],FUN="mean"),col="blue")
# lines(aggregate(pres~year,dat[dat$year<=2020,],FUN="min"),col="red")  #JS***
# plot(aggregate(abundance~year,dat,FUN="sum"),type="l",  lwd=2,ylab="Abundance", col="dark grey")
# lines(aggregate(abundance~year,dat[dat$year<=2020,],FUN="sum"),col="blue")
# plot(aggregate(temp~year,dat,FUN="min"),type="l",ylab="Temperature",ylim=c(-2,15), col="dark grey")
# lines(aggregate(temp~year,dat,FUN="max"),type="l",col="dark grey")
# lines(aggregate(temp~year,dat,FUN="mean"),type="l")

## JS: 
#randomly subset the data - simulates ~imperfect data collection
num_obs <- 500
dat_histc <- dat_hist[sample(1:nrow(dat_hist), num_obs, replace=F),]


#----Build GAM Models----
dat_histc$log_abundance <- log(dat_histc$abundance)
gam1.p <- gam(pres ~ s(temp,bs='gp') , data=dat_histc, family=binomial)
gam1.a <- gam(log_abundance ~ s(temp,bs='gp')  , data=dat_histc[dat_histc$abundance>0,], family=gaussian)

# saveRDS(gam1.p, paste0(Sim1,'GAM_Sim1_binom.rds'))
# saveRDS(gam1.a, paste0(Sim1,'GAM_Sim1_lognorm.rds'))
# # gam1.p <- readRDS( paste0(Sim1,'GAM_Sim1_binom.rds')) #read in model object if required
# # gam1.a <- readRDS( paste0(Sim1,'GAM_Sim1_lognorm.rds'))

# summary(gam1.p)
# summary(gam1.a)
# plot(gam1.p)
# plot(gam1.a)

##JS: Additional GAMs
# gam with default smoothness
M1.1 <- gam(round(abundance) ~ s(temp), data=dat_histc, family=poisson)
#summary(M1.1); plot(M1.1)

# gam with restricted smoothness
M1.2 <- gam(round(abundance) ~ s(temp, k=4), data=dat_histc, family=poisson)

# gam with restricted smoothness and zeros added at upper thermal limit (~8C)
dat_upper <- dat_histc[1:(nrow(dat_histc)*0.05),]  #add 5% extra rows as zeros  ***need a smart way to calculate penalty here; even very few data points can have big impact
dat_upper[] <- 0
dat_upper$temp <- 8  #estimated upper thermal limit
dat_upper$abundance <- 0  #all zeros
dat_histc2 <- rbind(dat_histc, dat_upper)
M1.3 <- gam(round(abundance) ~ s(temp, k=4), data=dat_histc2, family=poisson)

##JS: PLOT responses
par(mfrow=c(3,2))
ylim2 <- 70
new_dat <- data.frame(temp=seq(0,max(dat_hist$temp),length=100))
new_dat2 <- data.frame(temp=seq(0,7,length=100))
#actual TPC
xx <- seq(0, 7, length=100)
yy <- dnorm(xx, mean=3, sd=2)
plot(xx, yy, type="l", lty=2, main="Actual TPC", col="red", xlim=c(0,8), ylab="suitability", xlab="Temp", ylim=c(0,0.25))
xlim <- round(100*(max(dat_hist$temp)/7))
lines(xx[1:xlim], yy[1:xlim], lwd=2)
#delta model
P1p <- predict(gam1.p,new_dat2,type='response')
P1a <- predict(gam1.a,new_dat2,type="response")
P1 <- P1p*exp(P1a)
plot(new_dat2$temp, P1, type="l", main="delta-gam", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2),
     ylab="Abundance", xlab="Temp")
P2p <- predict(gam1.p,new_dat,type='response')
P2a <- predict(gam1.a,new_dat,type="response")
P2 <- P2p*exp(P2a)
points(dat_histc$temp, dat_histc$abundance, col="grey")
lines(new_dat$temp, P2, lwd=2)
#gam 1
plot(new_dat2$temp, predict(M1.1, newdata=new_dat2, type="response"), type="l",
     main="Poisson", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
points(dat_histc$temp, dat_histc$abundance, col="grey")
lines(new_dat$temp, predict(M1.1, newdata=new_dat, type="response"), lwd=2)
#gam 2
plot(new_dat2$temp, predict(M1.2, newdata=new_dat2, type="response"), type="l",
     main="Poisson, k=4", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
points(dat_histc$temp, dat_histc$abundance, col="grey")
lines(new_dat$temp, predict(M1.2, newdata=new_dat, type="response"), lwd=2)
#gam 3
plot(new_dat2$temp, predict(M1.3, newdata=new_dat2, type="response"), type="l",
     main="Poisson, k=4, upperTL", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
points(dat_histc2$temp, dat_histc2$abundance, col="grey")
lines(new_dat$temp, predict(M1.3, newdata=new_dat, type="response"), lwd=2)
par(mfrow=c(1,1))

#----Boosted Regression Tree----
#Make sure >1000 trees fitted
brt1.a <- gbm.step(data=dat_histc[dat_histc$abundance>0,], gbm.x = c("temp"),gbm.y = 'log_abundance',family = "gaussian",tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6)
brt1.p <- gbm.step(data=dat_histc, gbm.x = c("temp"),gbm.y = 'pres',family = "bernoulli",tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.6)
# saveRDS(brt1.a,paste0(Sim1,'BRT_Sim1_lognorm.rds'))
# saveRDS(brt1.p,paste0(Sim1,'BRT_Sim1_binom.rds'))
# # brt1.a <- readRDS(paste0(Sim1,'BRT_Sim1_lognorm.rds'))#read in model object if required
# # brt1.p <- readRDS(paste0(Sim1,'BRT_Sim1_binom.rds'))

# dev_eval=function(model_object){
#   null <- model_object$self.statistics$mean.null
#   res <- model_object$self.statistics$mean.resid
#   dev=((null - res)/null)*100 
#   return(dev)
# }
# dev_eval(brt1.p)
# dev_eval(brt1.a)
# 
# plot(brt1.p)
# plot(brt1.a)

#----Make Predictions for the future----
#GAM Hindcast (aka Fitted values)
dat_hist$gam1.p <- predict(gam1.p,dat_hist,type='response')
dat_hist$gam1.a <- predict(gam1.a,dat_hist,type="response")
dat_hist$gam1 <- dat_hist$gam1.p*exp(dat_hist$gam1.a)
dat_hist$brt1.p <- predict(brt1.p,dat_hist,n.trees=brt1.p$gbm.call$best.trees,type='response')
dat_hist$brt1.a <- predict(brt1.a,dat_hist,n.trees=brt1.a$gbm.call$best.trees,type='response')
dat_hist$brt1 <- dat_hist$brt1.p*exp(dat_hist$brt1.a)
dat_hist$gamP <- predict(M1.2, dat_hist, type="response")  #JS
dat_hist$gamPT <- predict(M1.3, dat_hist, type="response")  #JS

#GAM Forecast
dat_fcast$gam1.p <- predict(gam1.p,dat_fcast,type='response')
dat_fcast$gam1.a <- predict(gam1.a,dat_fcast,type="response")
dat_fcast$gam1 <- dat_fcast$gam1.p*exp(dat_fcast$gam1.a)
dat_fcast$brt1.p <- predict(brt1.p,dat_fcast,n.trees=brt1.p$gbm.call$best.trees,type='response')
dat_fcast$brt1.a <- predict(brt1.a,dat_fcast,n.trees=brt1.a$gbm.call$best.trees,type='response')
dat_fcast$brt1 <- dat_fcast$brt1.p*exp(dat_fcast$brt1.a)
dat_fcast$gamP <- predict(M1.2, dat_fcast, type="response")  #JS
dat_fcast$gamPT <- predict(M1.3, dat_fcast, type="response")  #JS

# #Standard errors Historical
# testCI1.a <- predict(gam1.a, dat_hist, type="response",se.fit=TRUE)
# testCI1.p <- predict(gam1.p, dat_hist, type="response",se.fit=TRUE)
# dat_hist$gam1.high <- exp((testCI1.a$fit + (testCI1.a$se.fit))) * (testCI1.p$fit + (testCI1.p$se.fit))
# dat_hist$gam1.low <- exp((testCI1.a$fit - (testCI1.a$se.fit))) * (testCI1.p$fit - (testCI1.p$se.fit))
# 
# #Standard errors Forecast
# testCI1.a <- predict(gam1.a, dat_fcast, type="response",se.fit=TRUE)
# testCI1.p <- predict(gam1.p, dat_fcast, type="response",se.fit=TRUE)
# dat_fcast$gam1.high <- exp((testCI1.a$fit + (testCI1.a$se.fit))) * (testCI1.p$fit + (testCI1.p$se.fit))
# dat_fcast$gam1.low <- exp((testCI1.a$fit - (testCI1.a$se.fit))) * (testCI1.p$fit - (testCI1.p$se.fit))

#Errors from BRT can be generated, I just haven't added the code yet (steph)

#----Compare abundance predictons----
#Quick and dirty plots (will convert to ggplot at some point)
par(mfrow=c(1,2))
#Historical patterns
plot(aggregate(abundance~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance", ylim=c(0,20000))
lines(aggregate(gam1~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="blue")
#lines(aggregate(gam1.high~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="light blue", lty=2)
#lines(aggregate(gam1.low~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="light blue", lty=2)
lines(aggregate(brt1~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="red")
lines(aggregate(gamP ~ year, dat_hist, FUN="sum"), lwd=2, col="green")  #JS
lines(aggregate(gamPT ~ year, dat_hist, FUN="sum"), lwd=2, col="purple")  #JS
legend("bottomleft",c("Truth","GAM-delta","BRT", "GAM-P", "GAM-PT"),lty=1,
       col=c("black","blue","red", "green", "purple"),bty="n")

#Future patterns
plot(aggregate(abundance~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance", ylim=c(0,20000))
lines(aggregate(gam1~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='blue')
#lines(aggregate(gam1.high~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='light blue')
#lines(aggregate(gam1.low~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='light blue')
lines(aggregate(brt1~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='red')
lines(aggregate(gamP ~ year, dat_fcast, FUN="sum"), lwd=2, col="green")  #JS
lines(aggregate(gamPT ~ year, dat_fcast, FUN="sum"), lwd=2, col="purple")  #JS
legend("bottomleft",c("Truth","GAM-delta","BRT", "GAM-P", "GAM-PT"),lty=1,
       col=c("black","blue","red", "green", "purple"),bty="n")



# #Historical COG
# cog_hist_lat <- as.data.frame(matrix(NA,nrow=20,ncol=6))
# colnames(cog_hist_lat) <- c("year","truth","gam1","gam1.p","brt1","brt1.p")
# counter=1
# for (y in 2001:2020){
#   cog_hist_lat[counter,1] <- y
#   cog_hist_lat[counter,2] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$abundance[dat_hist$year==y])
#   cog_hist_lat[counter,3] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$gam1[dat_hist$year==y])
#   cog_hist_lat[counter,4] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$gam1.p[dat_hist$year==y])
#   cog_hist_lat[counter,5] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$brt1[dat_hist$year==y])
#   cog_hist_lat[counter,6] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$brt1.p[dat_hist$year==y])
#   counter = counter + 1
# }
# head(cog_hist_lat)
# plot(cog_hist_lat$year,cog_hist_lat$truth, type='b')
# lines(cog_hist_lat$year,cog_hist_lat$gam1, type='b', col="blue")
# lines(cog_hist_lat$year,cog_hist_lat$brt1, type='b', col="red")
# 
# #Future COG
# cog_fcast_lat <- as.data.frame(matrix(NA,nrow=80,ncol=6))
# colnames(cog_fcast_lat) <- c("year","truth","gam1","gam1.p","brt1","brt1.p")
# counter=1
# for (y in 2021:2100){
#   cog_fcast_lat[counter,1] <- y
#   cog_fcast_lat[counter,2] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$abundance[dat_fcast$year==y])
#   cog_fcast_lat[counter,3] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam1[dat_fcast$year==y])
#   cog_fcast_lat[counter,4] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam1.p[dat_fcast$year==y])
#   cog_fcast_lat[counter,5] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$brt1[dat_fcast$year==y])
#   cog_fcast_lat[counter,6] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$brt1.p[dat_fcast$year==y])
#   counter = counter + 1
# }
# head(cog_fcast_lat)
# plot(cog_fcast_lat$year,cog_fcast_lat$truth, type='b')
# lines(cog_fcast_lat$year,cog_fcast_lat$gam1, type='b', col="blue")
# lines(cog_fcast_lat$year,cog_fcast_lat$brt1, type='b', col="red")
# 
# #-----Plot Surface Predictions-----
# #Future
# Y = 2021
# #Truth
# ggplot(dat_fcast[dat_fcast$year==Y,],aes(Lon,Lat))+
#   geom_tile(aes(fill=abundance)) +
#   theme_classic() +
#   ggtitle("Truth")+
#   labs(y="Latitude") +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous( expand = c(0, 0)) +
#   theme(legend.title=element_blank(),
#         plot.title = element_text(hjust=0.5),
#         # axis.text.x=element_blank(),
#         # axis.ticks=element_blank(),
#         # axis.title.x=element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)) +
#   scale_fill_viridis()
# 
# #Gam
# ggplot(dat_fcast[dat_fcast$year==Y,],aes(Lon,Lat))+
#   geom_tile(aes(fill=gam1)) +
#   theme_classic() +
#   ggtitle("GAM")+
#   labs(y="Latitude") +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous( expand = c(0, 0)) +
#   theme(legend.title=element_blank(),
#         plot.title = element_text(hjust=0.5),
#         # axis.text.x=element_blank(),
#         # axis.ticks=element_blank(),
#         # axis.title.x=element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)) +
#   scale_fill_viridis()
# 
# #BRT
# ggplot(dat_fcast[dat_fcast$year==Y,],aes(Lon,Lat))+
#   geom_tile(aes(fill=brt1)) +
#   theme_classic() +
#   ggtitle("BRT")+
#   labs(y="Latitude") +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous( expand = c(0, 0)) +
#   theme(legend.title=element_blank(),
#         plot.title = element_text(hjust=0.5),
#         # axis.text.x=element_blank(),
#         # axis.ticks=element_blank(),
#         # axis.title.x=element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)) +
#   scale_fill_viridis()
