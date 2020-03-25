
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
abund_enviro <- "poisson" #can be "lnorm_low" (SB); "lnorm_high" (EW); or "poisson" (JS)
PA_shape <- "logistic_prev" #can be "logistic" (SB); "logistic_prev","linear" (JS)
temp_spatial <- "matern" #can be "simple" (SB); or "matern" (EW)
temp_diff <- c(1,4,3,7) #specifies min and max temps at year 1 and year 100 (e.g. temp_diff=c(1,3,5,7) means year 1 varies from 1-3C and year 100 from 5-7C). For non-ROMS data. 
#dir <- "~/Dropbox/WRAP Location^3/Rasters_2d_monthly/" #directory where ROMS data is stored (on dropbox, email steph for access)

#Run this function
# dat <- SimulateWorld_ROMS(PA_shape = PA_shape, abund_enviro = abund_enviro, dir = dir ) #takes a few mins
#OR this function
dat <- SimulateWorld(temp_diff = temp_diff,  temp_spatial = temp_spatial, PA_shape = PA_shape, abund_enviro = abund_enviro) #takes a few minutes
saveRDS(dat, "dat.rds")

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


#----Build GAM Models----

gam1 <- gam(abundance ~ s(temp), data=dat_hist, family=poisson)
plot(gam1)


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
))

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

save(glmm1,glmm2,glmm3,glmm4,glmm5,glmm6, file = "glmms.rdata")
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
  pred_summary = dplyr::mutate(pred, lon_cell = ceiling(Lon/4),
                               lat_cell = ceiling(Lat/4)) %>% 
    dplyr::filter(year <= max_year) %>%
    dplyr::group_by(lon_cell, lat_cell, year) %>% 
    dplyr::summarize(mean_obs = mean(abundance),
                     mean_pred = sum(exp(est)))
  
  return(pred_summary)
} 

pred_summary <- predict_glmm(model = glmm1, max_year = 2026)

ggplot(pred_summary,aes(mean_pred,mean_obs)) + geom_point(alpha = 0.4) + 
  geom_abline(intercept = 0, slope = 1) + facet_wrap(~year)