
setwd("C:/Users/James.Thorson/Desktop/Git/WRAP_Location_CaseStudy")

#----Load Library & Function----
library(dismo)
library(gbm)
library(mgcv)
library(ggplot2)
library(viridis)
library(BBmisc)
library(neuralnet)

source('SimulatedWorld_Function.R') #load simulation function
#source('SimulatedWorld_ROMS_Function.R') #load ROMS simulation function

#-----Simulate data----

#Set parameters for functions
abund_enviro <- "lnorm_low" #can be "lnorm_low" (SB); "lnorm_high" (EW); or "poisson" (JS)
PA_shape <- "logistic_prev" #can be "logistic" (SB); "logistic_prev","linear" (JS)
temp_spatial <- "matern" #can be "simple" (SB); or "matern" (EW)
temp_diff <- c(1,4,3,7) #specifies min and max temps at year 1 and year 100 (e.g. temp_diff=c(1,3,5,7) means year 1 varies from 1-3C and year 100 from 5-7C). For non-ROMS data. 
#dir <- "~/Dropbox/WRAP Location^3/Rasters_2d_monthly/" #directory where ROMS data is stored (on dropbox, email steph for access)

#Run this function
# dat <- SimulateWorld_ROMS(PA_shape = PA_shape, abund_enviro = abund_enviro, dir = dir ) #takes a few mins
#OR this function
dat <- SimulateWorld(temp_diff = temp_diff,  temp_spatial = temp_spatial, PA_shape = PA_shape, abund_enviro = abund_enviro) #takes a few minutes

#make headers consistent (Steph needs to update functions to fix this)
colnames(dat)[1:2] <- c("Lon","Lat")


#Make some quick plots to explore the data
#All Years
par(mfrow=c(2,2))
plot(aggregate(suitability~year,dat,FUN="mean"),type="l", lwd=2, ylab="Suitability",col="dark grey")
lines(aggregate(suitability~year,dat[dat$year<=2020,],FUN="mean"),col="blue")
plot(aggregate(pres~year,dat,FUN="mean"),type="l", lwd=2,ylab="Presence",col="dark grey")
lines(aggregate(pres~year,dat[dat$year<=2020,],FUN="mean"),col="blue")
plot(aggregate(abundance~year,dat,FUN="sum"),type="l",  lwd=2,ylab="Abundance", col="dark grey")
lines(aggregate(abundance~year,dat[dat$year<=2020,],FUN="sum"),col="blue")
plot(aggregate(temp~year,dat,FUN="min"),type="l",ylab="Temperature",ylim=c(-2,15), col="dark grey")
lines(aggregate(temp~year,dat,FUN="max"),type="l",col="dark grey")
lines(aggregate(temp~year,dat,FUN="mean"),type="l")
par(mfrow=c(1,1))

plot(aggregate(temp~year,dat,FUN="mean"),type="l",ylab="Temperature", col="dark grey")


#Create dataframe with historical/forecast data
dat_hist <- dat[dat$year<=2020,]
dat_fcast <- dat[dat$year>2020,]


#---- GAMs ----

dat_hist$fYear <- as.factor(dat_hist$year)

gam_enviro <- gam(abundance ~ s(temp), data=dat_hist, family=tw())
plot(gam_enviro); AIC(gam_enviro)  #AIC=25340

gam_ST1 <- gam(abundance ~ s(temp) + s(Lat,Lon), data=dat_hist, family=tw())
plot(gam_ST1, pages=1); AIC(gam_ST1)  #25320

# gam_ST2 <- gam(abundance ~ s(temp) + s(Lat,Lon, by=fYear), data=dat_hist, family=tw())
# plot(gam_ST2)

gam_ST3 <- gam(abundance ~ s(temp) + te(Lat,Lon,year), data=dat_hist, family=poisson)
plot(gam_ST3)

gam_ST4 <- gamm(abundance ~ s(temp) + s(Lat,Lon), correlation=corGaus(.1, form=~Lat+Lon), data=dat_hist, family=poisson)
plot(gam_ST4$gam)


#---- GAMs with Metabolic 'COnstraint' ----

#add fake O2 and metabolic index to data
dat_hist$O2 <- rnorm(nrow(dat_hist), 3, 1)  #fake oxygen data
#dat$MI <- dat$suitability * dat$O2  #fake metabolic index data
dat_hist$MI <- dat_hist$suitability * exp(dat_hist$O2)  #fake metabolic index data


# fit gam with Gaussian process smoothers so variances are additive in log-space
gam1 <- gam(abundance ~ s(temp,bs='gp'), data=dat_hist, family=poisson(link=log))

#fit gam with metabolic index as offset; use log(MI) so that its multiplicative in natural space
gam2 <- gam(abundance ~ s(temp,bs='gp'), data=dat_hist, family=poisson(link=log), offset=log(dat_hist$MI) )

#fit gam with metabolic index as linear covariate (adds one parameter relative to gam2)
gam3 <- gam(abundance ~ s(temp,bs='gp') + dat_hist$MI, data=dat_hist, family=poisson(link=log) )

#fit gam with metabolic index as GP smoothed response (adds effective_degrees_of_freedome relative to gam2)
gam4 <- gam(abundance ~ s(temp,bs='gp') + s(dat_hist$MI,bs='gp'), data=dat_hist, family=poisson(link=log) )

#fit with spatially varying impact of linear metabolic inde
gam5 <- gam(abundance ~ s(temp,bs='gp') + s(Lon,Lat,by=dat_hist$MI), data=dat_hist, family=poisson(link=log) )

# see Degrees of freedom
anova(gam1,gam2,gam3,gam4,gam5) #, text="Chisq")


#---- GLMMs ----