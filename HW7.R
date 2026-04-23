#Jamie Boyd
#HW7

# read in greenhouse gas data from reservoirs
ghg <- read.csv("/cloud/project/activity07/Deemer_GHG_Data.csv")

install.packages(c("ggplot2", "dplyr", "olsrr", "PerformanceAnalytics", "lubridate", "forecast" ))

library(dplyr)
library(ggplot2)
library(olsrr)
library(PerformanceAnalytics)
library(lubridate)
library(forecast)

#In Class Activity/Questions ----

# log transform methane fluxes
ghg$log.ch4 <- log(ghg$ch4+1)
ghg$log.age <- log(ghg$age)
ghg$log.DIP <- log(ghg$DIP+1)
ghg$log.precip <- log(ghg$precipitation)

unique(ghg$Region)
# binary variable for boreal region
ghg$BorealV <- ifelse(ghg$Region == "Boreal",1,0)
# binary variable for tropical region
ghg$TropicalV <- ifelse(ghg$Region == "Tropical",1,0)


# binary variable for alpine region
ghg$AlpineV <- ifelse(ghg$Alpine == "yes",1,0)
# binary variable for known hydropower
ghg$HydroV <- ifelse(ghg$hydropower == "yes",1,0)

# multiple regression
# creates a model object
mod.full <- lm(log.ch4 ~ HydroV + airTemp+
                 log.age+mean.depth+
                 log.DIP+
                 log.precip+ BorealV, data=ghg) #uses the data argument to specify dataframe
summary(mod.full)

res.full <- rstandard(mod.full)
fit.full <- fitted.values(mod.full)

# qq plot
qqnorm(res.full, pch=19, col="grey50")
qqline(res.full)

# shapiro-wilks test
shapiro.test(res.full)
plot(fit.full,res.full, pch=19, col="grey50")
abline(h=0)

# isolate continuous model variables into data frame:

reg.data <- data.frame(ghg$HydroV, ghg$airTemp,
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
plot(full.step )

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

#Almond Orchard Tutorial ----

ETdat <- read.csv("/cloud/project/activity07/ETdata.csv")

unique(ETdat$crop)

# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(almond, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)

almondTrend <- almond_dec$trend
almondSeason <- almond_dec$seasonal

acf(na.omit(almond_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

pacf.plot <- pacf(na.omit(almond_ts))

almond_y <- na.omit(almond_ts)
model1 <- arima(almond_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model1

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
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

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

#Start of HW7 ----

#Question 1 ----

ghg$co2_trans = 1/(ghg$co2 + 1000)

#Create new model
mod.full.trans = lm(co2_trans ~ airTemp+
                      log.age+mean.depth+
                      ghg$chlorophyll.a+ ghg$HydroV+
                      log.DIP+
                      log.precip+ BorealV, data=ghg) #uses the data argument to specify dataframe
summary(mod.full.trans)
nobs(mod.full.trans)

#Check main assumptions
res.full_trans = rstandard(mod.full.trans)
fit.full_trans = fitted.values(mod.full.trans)

#Check normality of residuals
# qq plot
qqnorm(res.full_trans, pch=19, col="grey50")
qqline(res.full_trans)

# shapiro-wilks test
shapiro.test(res.full_trans)
plot(fit.full_trans,res.full_trans, pch=19, col="grey50")
abline(h=0)

#Checking for Multicolinearity

reg.data_trans <- data.frame(ghg$airTemp,
                        ghg$log.age,ghg$mean.depth,
                        ghg$chlorophyll.a,
                        ghg$log.DIP,
                        ghg$log.precip+ghg$HydroV)
chart.Correlation(reg.data, histogram=TRUE, pch=19) # make a correlation matrix


#Question 3 ----

# Almonds
# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

ggplot(almond, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)", title="Almonds")

almond_ts <- ts(almond$ET.in, start = c(2016,1), frequency=12)

almond_dec <- decompose(almond_ts)
plot(almond_dec)


# Pistachios
pistachio <- ETdat %>%
  filter(crop == "Pistachios") %>%
  group_by(date) %>%
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE))

ggplot(pistachio, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthly evapotranspiration (in)", title="Pistachios")

pistachio_ts <- ts(pistachio$ET.in, start = c(2016,1), frequency=12)
pistachio_dec <- decompose(pistachio_ts)
plot(pistachio_dec)


#Fallow
fallow <- ETdat %>%
  filter(crop == "Fallow/Idle Cropland") %>%
  group_by(date) %>%
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE))

ggplot(fallow, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthly evapotranspiration (in)", title="Fallow")

fallow_ts <- ts(fallow$ET.in, start = c(2016,1), frequency=12)
fallow_dec <- decompose(fallow_ts)
plot(fallow_dec)


#Corn
corn <- ETdat %>%
  filter(crop == "Corn") %>%
  group_by(date) %>%
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE))

ggplot(corn, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthly evapotranspiration (in)", title="Corn")

corn_ts <- ts(corn$ET.in, start = c(2016,1), frequency=12)
corn_dec <- decompose(corn_ts)
plot(corn_dec)


#Grapes, Grapes (Table/Raisin)
grapes <- ETdat %>%
  filter(crop == "Grapes (Table/Raisin)") %>%
  group_by(date) %>%
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE))

ggplot(grapes, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthly evapotranspiration (in)", title="Grapes")

grapes_ts <- ts(grapes$ET.in, start = c(2016,1), frequency=12)

grapes_dec <- decompose(grapes_ts)
plot(grapes_dec)

#Question 4 ----

#pistachios
acf(na.omit(pistachio_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

pacf.plot <- pacf(na.omit(pistachio_ts))

pistachio_y <- na.omit(pistachio_ts)
model1_p <- arima(pistachio_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model1_p

model4_p <- arima(pistachio_y, # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4_p

# calculate fit
AR_fit1_p <- pistachio_y - residuals(model1_p) 
AR_fit4_p <- pistachio_y - residuals(model4_p)
#plot data
plot(pistachio_y)
# plot fit
points(AR_fit1_p, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fit4_p, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

newP <- forecast(model4_p)
newP

#make dataframe for plotting
newPf <- data.frame(newP)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newPf$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = pistachio, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(pistachio$date[1]),newPf$dateF[24])+  # Plotting original data
  geom_line(data = newPf, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newPf, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)", title="Pistachio")


#fallow
#pistachios
acf(na.omit(fallow_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

pacf.plot <- pacf(na.omit(fallow_ts))

fallow_y <- na.omit(fallow_ts)
model1_f <- arima(fallow_y , # data 
                  order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model1_f

model4_f <- arima(fallow_y, # data 
                  order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4_f

# calculate fit
AR_fit1_f <- fallow_y - residuals(model1_p) 
AR_fit4_f <- fallow_y - residuals(model4_p)
#plot data
plot(fallow_y)
# plot fit
points(AR_fit1_f, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fit4_f, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

newf <- forecast(model4_f)
newf

#make dataframe for plotting
newff <- data.frame(newf)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newff$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = fallow, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(fallow$date[1]),newff$dateF[24])+  # Plotting original data
  geom_line(data = newff, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newff, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)", title="Fallow")
