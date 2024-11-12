install.packages(c("dplyr", "ggplot2", "olsrr", "PerformanceAnalytics"))
library(dplyr)
library(ggplot2)
library(olsrr)
library(PerformanceAnalytics)

# read in greenhouse gas data from reservoirs
ghg <- read.csv("/cloud/project/Deemer_GHG_Data.csv")

# log transform methane fluxes, age, dip, precip
ghg$log.ch4 <- log(ghg$ch4+1)
ghg$log.age <- log(ghg$age)
ghg$log.DIP <- log(ghg$DIP+1)
ghg$log.precip <- log(ghg$precipitation)

#make regions binary to investigate individual impacts
unique(ghg$Region) #check regions

# binary variable for boreal region
ghg$BorealV <- ifelse(ghg$Region == "Boreal",1,0)
# binary variable for tropical region
ghg$TropicalV <- ifelse(ghg$Region == "Tropical",1,0)
# binary variable for alpine
ghg$AlpineV <- ifelse(ghg$Alpine == "yes",1,0)
# binary variable for known hydropower
ghg$HydroV <- ifelse(ghg$hydropower == "yes",1,0)

# multiple regression
# creates a model object
mod.full <- lm(log.ch4 ~ airTemp+
                 log.age+mean.depth+
                 log.DIP+
                 log.precip+ BorealV, data=ghg) #uses the data argument to specify dataframe
summary(mod.full)

#assumption testing
res.full <- rstandard(mod.full)
fit.full <- fitted.values(mod.full)

#qq plot and qqline, normality of residuals
qqnorm(res.full, pch=19, col="grey50")
qqline(res.full)
#shapiro-wilks test, normality
shapiro.test(res.full)

#plot residuals
plot(fit.full,res.full, pch=19, col="grey50")
abline(h=0)

# isolate continuous model variables into data frame:
reg.data <- data.frame(ghg$airTemp,
                       ghg$log.age,ghg$mean.depth,
                       ghg$log.DIP,
                       ghg$log.precip)

# make a correlation matrix 
chart.Correlation(reg.data, histogram=TRUE, pch=19)

# run stepwise
full.step <- ols_step_forward_aic(mod.full)
# view table
full.step 
# check full model
full.step$model
# plot AIC over time
plot(full.step)

# prediction with interval for predicting a point
predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="prediction")

# look at prediction with 95% confidence interval of the mean
predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="confidence")

#in class prommpt 1
ETdat <- read.csv("/cloud/project/ETdata.csv")
unique(ETdat$crop)
install.packages(c("lubridate", "forecast"))
library(lubridate)
library(forecast)

# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm = T)) # average fields

# visualize the data
ggplot(almond, aes(x = ymd(date), y = ET.in)) +
  geom_point() +
  geom_line() +
  labs(x = "Year",  y= "Monthy Evapotranspiration (in)") +
  theme_classic()

# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)

# plot autocorrelations over lag period
acf(na.omit(almond_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)
# partial acf, can better represent potential lag periods
pacf.plot <- pacf(na.omit(almond_ts))

#run autoregressive modelwith more coefficients
almond_y <- na.omit(almond_ts)
model1 <- arima(almond_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model1

# change order
model4 <- arima(almond_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4

# calculate fit
AR_fit1 <- almond_y - residuals(model1) 
AR_fit4 <- almond_y - residuals(model4)
#plot data
plot(almond_y)
# plot fit
points(AR_fit1, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("right", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

# forecast data with uncertainty
newAlmond <- forecast(model4)
newAlmond

#make dataframe for plotting
newAlmondF <- data.frame(newAlmond)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newAlmondF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = almond, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(almond$date[1]),newAlmondF$dateF[24])+  # Plotting original data
  geom_line(data = newAlmondF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newAlmondF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")

#Question 1
#transform and design a regression analysis about the impact of reservoir characteristics on C fluxes
#consider the environmental conditions that impact carbon dioxide fluxes
#the availability of data, and the assumptions of ordinary least squares regression. 
#create a regression table including an R2 and the sample size 

#transform data
ghg$co2_transform <- 1/(ghg$co2 + 1000)

# multiple regression
# creates a model object
model.flux <- lm(co2_transform ~ airTemp +
                 log.age + mean.depth + surface.area + chlorophyll.a + log.DIP + 
                 log.ch4 + runoff + log.precip + BorealV + TropicalV, 
               data=ghg) #uses the data argument to specify dataframe
summary(model.flux)
#remove nas in new dataframe
model.flux.clean <- na.omit(model.flux)

#assumption testing
res.flux <- rstandard(model.flux.clean)
fit.flux <- fitted.values(model.flux.clean)

#qq plot and qqline, normality of residuals
qqnorm(res.flux, pch=19, col="grey50")
qqline(res.flux)
#shapiro-wilks test, normality
shapiro.test(res.flux)

#plot residuals
plot(fit.flux,res.flux, pch=19, col="grey50")
abline(h=0)

# isolate continuous model variables into data frame:
flux.data <- data.frame(ghg$airTemp, ghg$surface.area, ghg$chlorophyll.a,
                       ghg$log.age,ghg$mean.depth,
                       ghg$log.DIP, ghg$ch4, ghg$runoff,
                       ghg$log.precip)


# make a correlation matrix 
chart.Correlation(flux.data, histogram=TRUE, pch=19)

########## run stepwise
flux.step <- ols_step_forward_aic(model.flux.clean)
# view table
flux.step 
# check full model
flux.step$model
# plot AIC over time
plot(flux.step)

####### prediction with interval for predicting a point
predict.lm(model.flux.clean, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="prediction")

####### look at prediction with 95% confidence interval of the mean
predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="confidence")

#Question 2
#Decompose the evapotranspiration time series for almonds, pistachios, fallow/idle fields, corn, and table grapes. 
#Include plots of your decomposition.

#ALMOND
#average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm = T)) # average fields
# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit
# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)
# plot autocorrelations over lag period
acf(na.omit(almond_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)
almond.acf.plot <- acf(na.omit(almond_ts))

#PISTACHIOS
# average fields for each month for pistachio
pistachios <- ETdat %>% # ET data
  filter(crop == "Pistachios") %>% # only use pistachio fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm = T)) # average fields
# pistachio ET time series
pistachios_ts <- ts(pistachios$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit
# decompose pistachios ET time series
pistachios_dec <- decompose(pistachios_ts)
# plot decomposition
plot(pistachios_dec)
# plot autocorrelations over lag period
acf(na.omit(pistachios_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)
pistachios.acf.plot <- acf(na.omit(pistachios_ts))

#FALLOW/IDLE CROPLAND
# average fields for each month for fallow/idle
fallow <- ETdat %>% # ET data
  filter(crop == "Fallow/Idle Cropland") %>% # only use fallow/idle fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm = T)) # average fields
# fallow/idle ET time series
fallow_ts <- ts(fallow$ET.in, # data
                    start = c(2016,1), #start year 2016, month 1
                    #first number is unit of time and second is observations within a unit
                    frequency= 12) # frequency of observations in a unit
# decompose fallow/idle ET time series
fallow_dec <- decompose(fallow_ts)
# plot decomposition
plot(fallow_dec)
# plot autocorrelations over lag period
acf(na.omit(fallow_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)
fallow.acf.plot <- acf(na.omit(fallow_ts))

#CORN
# average fields for each month for corn
corn <- ETdat %>% # ET data
  filter(crop == "Corn") %>% # only use corn fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm = T)) # average fields
# corn ET time series
corn_ts <- ts(corn$ET.in, # data
                    start = c(2016,1), #start year 2016, month 1
                    #first number is unit of time and second is observations within a unit
                    frequency= 12) # frequency of observations in a unit
# decompose corn ET time series
corn_dec <- decompose(corn_ts)
# plot decomposition
plot(corn_dec)
# plot autocorrelations over lag period
acf(na.omit(corn_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)
corn.acf.plot <- acf(na.omit(corn_ts))

#GRAPES (TABLE/RAISIN)
# average fields for each month for grapes
grapes <- ETdat %>% # ET data
  filter(crop == "Grapes (Table/Raisin") %>% # only use grape fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm = T)) # average fields
# corn ET time series
grapes_ts <- ts(corn$ET.in, # data
              start = c(2016,1), #start year 2016, month 1
              #first number is unit of time and second is observations within a unit
              frequency= 12) # frequency of observations in a unit
# decompose grapes ET time series
grapes_dec <- decompose(grapes_ts)
# plot decomposition
plot(grapes_dec)
# plot autocorrelations over lag period
acf(na.omit(grapes_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)
grapes.acf.plot <- acf(na.omit(grapes_ts))

#Question 3
#Design an autoregressive model for pistachios and fallow/idle fields. 
#Forecast future evapotranspiration for each field
#Make a plot that includes historical and forecasted evapotranspiration for the crops

#PISTACHIOS
#run autoregressive modelwith more coefficients
pistachios_y <- na.omit(pistachios_ts)
model4.pist <- arima(pistachios_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4.pist
# calculate fit
AR_fit4_pist <- pistachios_y - residuals(model4.pist)
#plot data
plot(pistachios_y)
# plot fit
points(AR_fit4_pist, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
# forecast data with uncertainty
newPistachios <- forecast(model4.pist)
newPistachios
#make dataframe for plotting
newPistachiosF <- data.frame(newPistachios)
# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newPistachiosF$dateF <- ymd(paste(years,"/",month,"/",1))
# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = pistachios, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(pistachios$date[1]),newPistachiosF$dateF[24])+  # Plotting original data
  geom_line(data = newPistachiosF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newPistachiosF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(title="Pistachios", x="Year", y="Evapotranspiration (in)")

#FALLOW/IDLE
#run autoregressive model with more coefficients
fallow_y <- na.omit(fallow_ts)
model4.fallow <- arima(fallow_y , # data 
                     order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4.fallow
# calculate fit
AR_fit4_fallow <- fallow_y - residuals(model4.fallow)
#plot data
plot(fallow_y)
# plot fit
points(AR_fit4_fallow, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
# forecast data with uncertainty
newFallow <- forecast(model4.fallow)
newFallow
#make dataframe for plotting
newFallowF <- data.frame(newFallow)
# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newFallowF$dateF <- ymd(paste(years,"/",month,"/",1))
# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = fallow, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(fallow$date[1]),newFallowF$dateF[24])+  # Plotting original data
  geom_line(data = newFallowF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newFallowF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(title="Fallow/Idle Cropland", x="Year", y="Evapotranspiration (in)")
