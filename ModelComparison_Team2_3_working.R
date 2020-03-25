
#setwd("C:/Users/James.Thorson/Desktop/Git/WRAP_Location_CaseStudy")

#----Load Library & Function----
library(dismo)
library(gbm)
library(mgcv)
library(ggplot2)
library(viridis)
library(BBmisc)
library(neuralnet)
devtools::install_github("pbs-assess/sdmTMB")
library(sdmTMB)
library(sp)
library(dplyr)

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
#saveRDS(dat, "dat.rds")

#make headers consistent (Steph needs to update functions to fix this)
colnames(dat)[1:2] <- c("Lon","Lat")

# Prepare to project coordinates (if using actual west coast grid (ROMS))
dat_ll <- dat
coordinates(dat_ll) <- c("Lon", "Lat")
proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
# convert to UTM with spTransform
dat_utm <- spTransform(dat_ll, 
                       CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
dat = as.data.frame(dat_utm) # convert back from sp object to data frame


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

gam_enviro <- gam(abundance ~ s(temp), data=dat_hist, family=tw(link="log"), method="REML")
plot(gam_enviro); AIC(gam_enviro)  #AIC=25340, p=1.213

gam_ST1 <- gam(abundance ~ s(temp) + s(Lat,Lon), data=dat_hist, family=tw(link="log"), method="REML")
plot(gam_ST1, pages=1); AIC(gam_ST1)  #25320, p=1.213

# gam_ST2 <- gam(abundance ~ s(temp) + s(Lat,Lon, by=fYear), data=dat_hist, family=Tweedie(p=1.213, link="log"))
# plot(gam_ST2)

gam_ST3 <- gam(abundance ~ s(temp) + te(Lat,Lon,year), data=dat_hist, family=Tweedie(p=1.213, link="log"))
plot(gam_ST3); AIC(gam_ST3)  #25293

gam_ST4 <- gamm(abundance ~ s(temp), correlation=corGaus(form=~Lat+Lon|fYear), data=dat_hist,
                family=Tweedie(p=1.213, link="log"))  #takes a long time to fit...
plot(gam_ST4$gam)

# dat_upper <- dat_hist[1:(nrow(dat_hist)*0.05),]  #add 5% extra rows as zeros  ***need a smart way to calculate penalty here; even very few data points can have big impact
# dat_upper[] <- 0
# dat_upper$temp <- 8  #estimated upper thermal limit
# dat_upper$abundance <- 0  #all zeros
# dat_hist2 <- rbind(dat_hist, dat_upper)
# gam_UTL <- gam(abundance ~ s(temp), data=dat_hist2, family=tw(link="log"), method="REML")
# 
# 
# #---- GAMs with Metabolic 'Constraint' ----
# 
# #add fake O2 and metabolic index to data
# dat_hist$O2 <- rnorm(nrow(dat_hist), 3, 1)  #fake oxygen data
# #dat$MI <- dat$suitability * dat$O2  #fake metabolic index data
# dat_hist$MI <- dat_hist$suitability * exp(dat_hist$O2)  #fake metabolic index data
# #create fake temp-dependent aeorbic scope data
# dat_hist$AScope <- dnorm(dat_hist$temp, 4, 1)  #pretend that we know exactly how temp affects performance, and it's perfectly correlated with spatial association
# 
# 
# ## Try with MI, or some temp-related index (Aerobic Scope?)
# # fit gam with Gaussian process smoothers so variances are additive in log-space
# gam1 <- gam(abundance ~ s(temp,bs='gp'), data=dat_hist, family=tw(link=log))
# plot(gam1)
# 
# #fit gam with metabolic index as offset; use log(MI) so that its multiplicative in natural space
# gam2 <- gam(abundance ~ s(temp,bs='gp'), data=dat_hist, family=tw(link=log), offset=log(dat_hist$AScope))
# plot(gam2)
# 
# #fit gam with metabolic index as linear covariate (adds one parameter relative to gam2)
# gam3 <- gam(abundance ~ s(temp,bs='gp') + AScope, data=dat_hist, family=tw(link=log) )
# plot(gam3)
# 
# #fit gam with metabolic index as GP smoothed response (adds effective_degrees_of_freedome relative to gam2)
# gam4 <- gam(abundance ~ s(temp,bs='gp') + s(AScope,bs='gp'), data=dat_hist, family=tw(link=log) )
# 
# #fit with spatially varying impact of linear metabolic inde
# gam5 <- gam(abundance ~ s(temp,bs='gp') + s(Lon,Lat,by=MI), data=dat_hist, family=tw(link=log) )
# 
# # see Degrees of freedom
# anova(gam1,gam2,gam3,gam4,gam5) #, text="Chisq")
# 
# 
# # ---- Compare some future predictions of temperature ----
# ylim2 <- 10
# new_dat <- data.frame(temp=seq(0,max(dat_hist$temp),length=100))
# new_dat$AScope <- dnorm(new_dat$temp, 4, 1)
# new_dat2 <- data.frame(temp=seq(0,7,length=100))
# new_dat2$AScope <- dnorm(new_dat2$temp, 4, 1)
# 
# par(mfrow=c(2,1))
# #actual TPC
# xx <- seq(0, 7, length=100)
# yy <- dnorm(xx, mean=4, sd=1)  #Must match function in SimulatedWorld function
# plot(xx, yy, type="l", lty=2, main="Actual TPC", col="red", xlim=c(0,8), ylab="suitability", xlab="Temp")
# xlim <- round(100*(max(dat_hist$temp)/7))
# lines(xx[1:xlim], yy[1:xlim], lwd=2)
# #gam enviro
# plot(new_dat2$temp, predict(gam_enviro, newdata=new_dat2, type="response"), type="l",
#      main="Enviro GAM", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
# points(dat_hist$temp, dat_hist$abundance, col="grey")
# lines(new_dat$temp, predict(gam_enviro, newdata=new_dat, type="response"), lwd=2)
# #gam 2
# plot(new_dat2$temp, predict(gam2, newdata=new_dat2, type="response"), type="l",  
#      main="Enviro Gam w AS Offset", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
# points(dat_hist$temp, dat_hist$abundance, col="grey")
# lines(new_dat$temp, predict(gam2, newdata=new_dat, type="response"), lwd=2) 
# #gam 3
# plot(new_dat2$temp, predict(gam3, newdata=new_dat2, type="response"), type="l",
#      main="Enviro Gam w AS Linear", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
# points(dat_hist$temp, dat_hist$abundance, col="grey")
# lines(new_dat$temp, predict(gam3, newdata=new_dat, type="response"), lwd=2)
# #gam 4
# plot(new_dat2$temp, predict(gam4, newdata=new_dat2, type="response"), type="l",
#      main="Enviro Gam w AS Smoother", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
# points(dat_hist$temp, dat_hist$abundance, col="grey")
# lines(new_dat$temp, predict(gam4, newdata=new_dat, type="response"), lwd=2)
# #gam UTL
# plot(new_dat2$temp, predict(gam_UTL, newdata=new_dat2, type="response"), type="l",
#      main="Enviro Gam w UTL", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
# points(dat_hist2$temp, dat_hist2$abundance, col="grey")
# lines(new_dat$temp, predict(gam_UTL, newdata=new_dat, type="response"), lwd=2)


#----Build GLMM Models----

# note that there is a much cleaner way of looping through models 
# but for now simplifying for clarity and consistency with GAM fit formatting above

# make SPDE (stochastic partial differential equation that these INLA-based methods rely on)
spde <- try(make_spde(x = dat$Lon, y = dat$Lat, n_knots = 50), silent=TRUE) # increase knots to ~250 for WC data
plot_spde(spde) # can vary number of knots to modify the mesh until you get what you want

# need to remove future years from test data set via weights
weights = rep(0,nrow(dat))
weights[which(dat$year <= 2020)] = 1

# model without covariates but with spatiotemporal random fields as AR1
glmm1 <- try(sdmTMB(
  formula = abundance ~ -1,
  time_varying = NULL,
  spde = spde,
  time = "year",
  family = tweedie(link = "log"),
  data = dat,
  anisotropy = TRUE,
  ar1_fields = TRUE,
  weights = weights,
  spatial_only = FALSE,
  spatial_trend = FALSE
))  #takes 5+ mins

# model without covariates but with spatial trend random field, and spatiotemporal random fields as IID
glmm2 <- try(sdmTMB(
  formula = abundance ~ -1,
  time_varying = NULL,
  spde = spde,
  time = "year",
  family = tweedie(link = "log"),
  data = dat,
  anisotropy = TRUE,
  ar1_fields = FALSE,
  weights = weights,
  spatial_only = FALSE,
  spatial_trend = TRUE
  ))

# model with environmental covariate but no spatiotemporal random fields
glmm3 <- try(sdmTMB(
  formula = abundance ~ -1 + temp + I(temp^2),
  time_varying = NULL,
  spde = spde,
  time = "year",
  family = tweedie(link = "log"),
  data = dat,
  anisotropy = TRUE,
  ar1_fields = FALSE,
  weights = weights,
  spatial_only = TRUE,
  spatial_trend = FALSE,
  quadratic_roots = TRUE
))

# model with environmental covariate, spatial trend random field, but no spatiotemporal random fields
glmm4 <- try(sdmTMB(
  formula = abundance ~ -1 + temp + I(temp^2),
  time_varying = NULL,
  spde = spde,
  time = "year",
  family = tweedie(link = "log"),
  data = dat,
  anisotropy = TRUE,
  ar1_fields = FALSE,
  weights = weights,
  spatial_only = TRUE,
  spatial_trend = TRUE,
  quadratic_roots = TRUE
))

# model with environmental covariate and spatiotemporal random fields as AR1
glmm5 <- try(sdmTMB(
  formula = abundance ~ -1 + temp + I(temp^2),
  time_varying = NULL,
  spde = spde,
  time = "year",
  family = tweedie(link = "log"),
  data = dat,
  anisotropy = TRUE,
  ar1_fields = TRUE,
  weights = weights,
  spatial_only = FALSE,
  spatial_trend = FALSE,
  quadratic_roots = TRUE
))

# model with environmental covariate spatial trend random field, and spatiotemporal random fields as IID
glmm6 <- try(sdmTMB(
  formula = abundance ~ -1 + temp + I(temp^2),
  time_varying = NULL,
  spde = spde,
  time = "year",
  family = tweedie(link = "log"),
  data = dat,
  anisotropy = TRUE,
  ar1_fields = FALSE,
  weights = weights,
  spatial_only = FALSE,
  spatial_trend = TRUE,
  quadratic_roots = TRUE
))

# note glmm5 and glmm6 not converging with current operating model output


# make prediction-- have to add dummy years as temporary adjustment to sdmTMB requirement
# that all time steps in original model must also be present in prediction data frame
predict_glmm <- function(model, max_year){
  dummy = data.frame(year = unique(dat$year),
                     Lat=dat$Lat[1],
                     Lon=unique(dat$Lon[1]),
                     temp = dat$temp[1])
  pred = predict(model, 
                 newdata=rbind(dat_fcast[,c("year","Lon","Lat","temp")], 
                               dummy), xy_cols = c("Lon", "Lat"))
  pred = pred[-c(seq(nrow(pred)-nrow(dummy)+1,nrow(pred))),] # drop dummy data
  pred$abundance = dat_fcast$abundance
  
  # aggregate predictions and observations at some coarse spatial resolution
  # to get rid of occurrence of 0s (if they arise from operating model)
  pred_summary = dplyr::mutate(pred, lon_cell = ceiling(Lon/2),
                               lat_cell = ceiling(Lat/2)) %>% 
    dplyr::filter(year <= max_year) %>%
    dplyr::group_by(lon_cell, lat_cell, year) %>% 
    dplyr::summarize(mean_obs = mean(abundance),
                     mean_pred = sum(exp(est)))
  
  ggplot(pred_summary,aes(mean_pred,mean_obs)) + geom_point(alpha = 0.1) + 
    geom_abline(intercept = 0, slope = 1) + facet_wrap(~year)
} 

predict_glmm(model = glmm1, max_year = 2026)

# plot(dat_fcast$abundance, exp(pred$est), xlab="observed", ylab="predicted")
# abline(0,1, col='red')



# ---- Compare prediction with actual future distribution ----
RMSE = function(p, o){
  sqrt(mean((p - o)^2))
}

par(mfrow=c(4,1), mar=c(3,4,2.5,1))

dat_fcast$gam_enviro <- predict(gam_enviro,dat_fcast,type="response")
plot(aggregate(abundance~year,dat_fcast,FUN="sum"),type="l", lwd=2, ylab="Abundance", ylim=c(0,2500), 
     xlim=c(2021,2080), xlab="", main='Gamm-enviro, Future Actual vs Predicted')
lines(aggregate(gam_enviro~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='blue')

save_error <- data.frame(year=seq(2021,2080), cor=0, rmse=0)
for (yy in 2021:2080) {
  dat_fcast_x <- dat_fcast[dat_fcast$year==yy,]
  save_error$cor[save_error$year==yy] <- cor(dat_fcast_x$abundance, dat_fcast_x$gam_enviro)
  save_error$rmse[save_error$year==yy] <- RMSE(dat_fcast_x$abundance, dat_fcast_x$gam_enviro)
}
plot(save_error$year, save_error$cor, type="l", ylab="Correlation", xlab="",
     main="Correlation, Predicted and Actual across domain")
plot(save_error$year, save_error$rmse, type="l", ylab="Correlation", xlab="",
     main="RMSE, Predicted and Actual across domain")

## COG
cog_fcast_lat <- as.data.frame(matrix(NA,nrow=nrow(dat_fcast),ncol=3))
colnames(cog_fcast_lat) <- c("year","truth","gam_enviro")
counter=1
for (y in unique(dat_fcast$year)){
  cog_fcast_lat[counter,1] <- y
  cog_fcast_lat[counter,2] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$abundance[dat_fcast$year==y])
  cog_fcast_lat[counter,3] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam_enviro[dat_fcast$year==y])
  counter = counter + 1
}
head(cog_fcast_lat)
plot(cog_fcast_lat$year,cog_fcast_lat$truth, type='b', ylab="COG 'latitude'", xlim=c(2021,2080),
     main="COG")
lines(cog_fcast_lat$year,cog_fcast_lat$gam_enviro, type='b', col="blue")




