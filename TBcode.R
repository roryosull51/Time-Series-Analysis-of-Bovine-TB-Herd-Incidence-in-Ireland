###########################################################################
## Project: TSA Bovine TB
## Script purpose: R code for the project
## Date: 28/04/24
## Author: Rory O Sullivan 20232721
###########################################################################

library(TSA)
library(tidyverse)
library(tseries)
library(fpp2)
library(forecast)


df <- read_csv(("TBincidence.csv"))


head(df)
df <- df[, !names(df) %in% c("Statistic Label", "Regional Veterinary Offices", "UNIT")]
head(df)
tail(df)
dim(df)

training <- df[-c(48:52), ]
head(training)
tail(training)
# 52 entries

####################### # Changing to TS
# Convert dataframe to time series
incidence_ts <- ts(df$VALUE, frequency = 4, start = c(2011, 1))

training_ts <- ts(training$VALUE, frequency = 4, start = c(2011, 1))


# Print the time series
print(incidence_ts)
plot(incidence_ts, col = "#69b3a2", main = "TB Herd Incidence (%) Time Series",
     xlab = "Year", ylab = "Herd Incidence (%)")
points(y = incidence_ts, x = time(incidence_ts),
       pch = as.vector(season(incidence_ts)), col = "royalblue3", cex = 0.8)
lines(training_ts, col = "black")


##################################################################
# Working with the training data 
##################################################################


################### Breaking down the model
training_tsDecom <- decompose(training_ts)
plot(training_tsDecom)

diff(range(training_ts))
diff(range(training_tsDecom$trend,na.rm=T))
diff(range(training_tsDecom$seasonal,na.rm=T))
diff(range(training_tsDecom$random,na.rm=T))

# Seems to be a quadratic trend going on
#### try fit quadratic here 
ti <- as.vector(time(training_ts))
linear <- lm(training_ts ~ ti)
model <- lm(training_ts~ ti + I(ti^2))
plot(training_ts, pch = as.vector(season(training_ts)), col = 4, cex = 0.8,
     ylab = "Herd Incidence (%)", xlab = "Time")
points(y = training_ts, x = time(training_ts),
       pch = as.vector(season(training_ts)), col = 4, cex = 0.8)
lines(x=ti, y=predict(model), lty=2, col=2)
lines(x=ti, y=predict(linear), lty=2, col=1)
summary(model)

# Line seems to capture the trend well.

################################################################
# Looking at Making the series stationary
################################################################

BC <-BoxCox.ar(training_ts, lambda=seq(0,2,0.1) )
BC$mle
BC$ci

# BoxCox returns MLE = 0.6 with wide CI of (0,2)
training_tsBC <- training_ts^0.6
plot(training_tsBC)
plot(training_ts)

d.bc <- diff(training_tsBC, lag = 4)
plot(d.bc)
d1.bc <- diff(d.bc)
plot(d1.bc)
d2.bc <- diff(d1.bc)
plot(d2.bc)
# Does nothing here really


# try Log transform
training_tsLog <- log(training_ts)
plot(training_ts)
plot(training_tsLog)

# Similarly does little here, leaving series still stationary
plot(training_tsLog)

# Look at Differencing
training_tsDiff <- diff(training_ts)
plot(training_tsDiff,  main="Difference-TB", 
     xlab="Year", ylab="Difference-TB")
# Seems to have removed the general trend anyways

############## Get rid of seasonal effects now
d1.train <- diff(training_ts)
d1.train.lag <- diff(d1.train, lag = 4)
plot(d1.train.lag,  main="Once Difference-TB", 
     xlab="Year", ylab="Difference-TB")


d2.train <- diff(training_ts, differences = 2)
d2.train.lag <- diff(d2.train, lag = 4)
plot(d2.train.lag,  main="Twice Differenced-TB", 
     xlab="Year", ylab="Twice Differenced-TB")


########## Difference log
d1.train.lag.log <- diff(diff(log(training_ts), lag = 4))
plot(d1.train.lag.log,  main="Difference-Log-TB", 
     xlab="Year", ylab="Difference-TB")
# No need for logging



###########################################################
# Now looking at plotting some ARIMA Models
##########################################################

################################ Linear Model ##########################

plot(d1.train.lag,  main="Linear Trend Model", 
     xlab="Year", ylab="diff(ts_d1, lag = 4)")
acf(d1.train.lag, lag.max = 20)
# MA(1) seasonal and thats all maybe
pacf(d1.train.lag, lag.max = 20)
# Decaying so MA model here strictly perhaps

plot(acf(d1.train.lag), main = "ACF plot",
    ylab = "ACF")
plot(pacf(d1.train.lag), main = "Partial ACF plot",
    ylab = "Partial ACF")
adf.test(d1.train.lag)



eacf(d1.train.lag)
# Wall at MA(3) indicative of the seasonality
plot(armasubsets(d1.train.lag, nar = 8, nma = 8))
# Pointing out importance of an MA(3) and even MA(5) non seasonal here

# initially checking SARIMA (0,1,5), (0,1,1)4 
linmod <- Arima(training_ts, order = c(0,1,1),
                seasonal = list(order = c(0,1,1), period = 4))
linmod

# Insig. MA(2,3) non S and S (MA(1)) AIC 13.47
# Final Model gone for (0,1,1), (0,1,1) AIC = 8.05

residlin <- rstandard(linmod)
fitlin <- fitted(linmod)
plot(as.vector(fitlin), as.vector(residlin), xlab = "Fitted Values", 
     ylab = "Residuals")
qqnorm(residlin)
qqline(residlin)
shapiro.test(residlin)

plot(residlin)
ggtsdisplay(stat2)
checkresiduals(linmod)
  
plot(hist(residlin), main = "Histogram of Residuals", xlab = "Residual")


plot(acf(residlin), main = "Residual ACF")

runs(residlin)

LB.test(linmod, lag = 15)
LBpvals <- rep(NA, 15)
for(i in 3:15){
  LBpvals[i] <- LB.test(linmod, lag=i)$p.value
}
plot(LBpvals, ylim=c(0,1)); abline(h=0.05, lty=2)

predict(linmod, n.ahead = 5)


################################ Box Cox Linear ##########################

plot(d1.bc)
acf(d1.bc, lag.max = 20)
# MA(1) seasonal and thats all maybe
pacf(d1.bc, lag.max = 20)
# Decaying so MA model here strictly perhaps


eacf(d1.bc)
# Wall at MA(3) indicative of the seasonality
plot(armasubsets(d1.bc, nar = 8, nma = 8))
# Pointing out importance of an MA(3) and even MA(1) non seasonal here
# Not pointing out MA seasonal, is pointing out AR seasonal and non-s

# initially checking SARIMA (0,1,5), (0,1,1)4 
boxmod <- Arima(training_tsBC, order = c(0,1,1),
                seasonal = list(order = c(0,1,1), period = 4))
boxmod

# Insig. MA(2,3) non S and S (MA(1)) AIC 13.47
# Final Model gone for (0,1,1), (0,1,1) AIC = - 75.23

residbox <- rstandard(boxmod)
fitbox <- fitted(boxmod)
plot(as.vector(fitbox), as.vector(residbox))
qqnorm(residbox)
qqline(residbox)
shapiro.test(residbox)


acf(residbox)


LBpvals <- rep(NA, 15)
for(i in 5:15){
  LBpvals[i] <- LB.test(boxmod, lag=i)$p.value
}
plot(LBpvals, ylim=c(0,1)); abline(h=0.05, lty=2)

predict(boxmod^1/0.6, n.ahead = 5)

# Predicting Box model
predictions <- predict(boxmod, n.ahead = 5)
(predictions$pred)^1/0.6


################################ Box Cox Quadratic ##########################

plot(d2.bc)
acf(d2.bc, lag.max = 20)
# MA(1) seasonal and maybe MA(1) non s
pacf(d2.bc, lag.max = 20)
# Decaying so MA model here strictly perhaps
# Try AR(2) non s and maybe an SAR(1)


eacf(d2.bc)
# Wall at MA(3) indicative of the seasonality
plot(armasubsets(d2.bc, nar = 8, nma = 8))
# AR 1, 6 non S. MA 1 5,6 non S

# initially checking SARIMA (1,2,2), (1,1,1)4 
boxmodquad <- Arima(training_tsBC, order = c(0,2,2),
                seasonal = list(order = c(0,1,1), period = 4))
boxmodquad
# Insig. AR1, ma2 and sar1
# Final Model gone for (0,2,2), (0,1,1) AIC = - 57.66
residboxquad <- rstandard(boxmodquad)
fitboxquad <- fitted(boxmodquad)
plot(as.vector(fitboxquad), as.vector(residboxquad))
qqnorm(residboxquad)
qqline(residboxquad)
shapiro.test(residboxquad)


acf(residbox)

LBpvals <- rep(NA, 15)
for(i in 5:15){
  LBpvals[i] <- LB.test(boxmodquad, lag=i)$p.value
}
plot(LBpvals, ylim=c(0,1)); abline(h=0.05, lty=2)

# Predicting Box model
predictions <- predict(boxmodquad, n.ahead = 5)
(predictions$pred)^1/0.6

################################ Quadratic Model #######################

# Due to quadratic trend, we must difference the series twice
plot(d2.train.lag, main = "Time Series Plot",
     ylab = "diff(ts_d2, lag = 4)")
acf(d2.train.lag)
# MA(1) Seasonal and MA(1) non-seasonal
pacf(d2.train.lag)
# Tapering Seasonal Lags, maybe an AR(2) non-seasonal though maybe tapering

#plot(acf(d2.train.lag), main = "ACF plot",
#    ylab = "ACF")
#plot(pacf(d2.train.lag), main = "Partial ACF plot",
#    ylab = "Partial ACF")

adf.test(d2.train.lag)
qqnorm(d2.train.lag)
qqline(d2.train.lag)


plot(armasubsets(d2.train.lag, nar =8, nma = 8))
# Indicating a non-seasonal AR(1) and MA(1)
# it isnt picking up on the seasonal MA



# Thinking ARIMA (2,2,1), (0,1,1)
quadmod <- Arima(training_ts, order = c(1,2,2),
                    seasonal = list(order = c(0,1,1), period = 4))
quadmod
# AR(2) non-seasonal isnt sig.
# Removing yields AR(1) non seasonal also insig.

quadmod2 <- arima(training_ts, order = c(0,2,2),
                 seasonal = list(order = c(0,1,1), period = 4))
quadmod2


############################## Examining its residuals
resid <- rstandard(quadmod2)
fit <- fitted(quadmod2)
plot(as.vector(fit), as.vector(resid), xlab = "Fitted Values", 
     ylab = "Residuals")

qqnorm(resid)
qqline(resid)
shapiro.test(resid)

plot(hist(resid), main = "Histogram of Residuals")

plot(acf(resid), main = "Residual ACF")

LBpvals <- rep(NA, 15)
for(i in 5:15){
  LBpvals[i] <- LB.test(quadmod, lag=i)$p.value
}
plot(LBpvals, ylim=c(0,1)); abline(h=0.05, lty=2)


# Predicting Quad model
predict(quadmod, n.ahead = 5)

predict(quadmod2, n.ahead = 5)

plot(predict(quadmod2, n.ahead = 5))

plot(quadmod2, n.ahead =5)
plot(forecast(quadmod2, h =5))
plot(forecast(quadmod, h =5))

################################ Log Model #############################

plot(d1.train.lag.log)
acf(d1.train.lag.log, lag.max = 20)
# MA(1) seasonal and maybe MA(1) non-seasonal
pacf(d1.train.lag.log, lag.max = 20)
# Decaying so MA model here strictly perhaps non-seasonal AR(1)


eacf(d1.train.lag.log)
# Wall at MA(3) indicative of the seasonality
plot(armasubsets(d1.train.lag.log, nar = 8, nma = 8))
# Pointing out importance of an MA(3) non seasonal here

checkresiduals(d1.train.lag.log)

# Thinking SARIMA (1,1,3), (0,1,1) 4
logmod <- Arima(log(training_ts), order = c(0,1,1),
                seasonal = list(order = c(0,1,1), period = 4))
logmod
# Non-Seasonal MA(2 and 3) insig
# this renders the non-seasonal AR(1) insig.
# Overall this model comes out with the lowest AIC

residlog <- rstandard(logmod)
fitlog <- fitted(logmod)
plot(as.vector(fitlog), as.vector(residlog))
qqnorm(residlog)
qqline(residlog)
shapiro.test(residlog)



acf(residlin)
LBpvals <- rep(NA, 15)
for(i in 5:15){
  LBpvals[i] <- LB.test(linmod, lag=i)$p.value
}
plot(LBpvals, ylim=c(0,1)); abline(h=0.05, lty=2)


# Predicting Log model
predictions <- predict(logmod, n.ahead = 5, transform = exp)
exp(predictions$pred)


######################################################
# Comparing the models on One Plot for Training Set
######################################################

plot(training_ts, col = "black", lty = 2, main = "TB Herd Incidence (%) Time Series",
     xlab = "Year", ylab = "Herd Incidence (%)", ylim = c(2,4.5))
#points(y = training_ts, x = time(training_ts),
#       pch = as.vector(season(training_ts)), col = "blue", cex = 0.8)
# Adding the good Quadratic Model
lines(quadmod2$fitted, col = "firebrick3")
#points(y = quadmod$fitted, x = time(quadmod$fitted),
#       pch = as.vector(season(incidence_ts)), col = "black", cex = 0.8)
# Adding the Linear Model(Linear and Logistic proved similar)
lines(linmod$fitted, col = "royalblue2")
# Adding log model?
lines(exp(logmod$fitted), col = "darkorange2")

######################################################
# Comparing the models on One Plot for Entire Set
######################################################

# First need to make new df for each model with the predicted values in them also

################### Linear Model
predict_lin <- predict(linmod, n.ahead = 5)
predict_lin$pred
linmod$fitt
lin.df <- data.frame(VALUE = c(linmod$fitted, predict_lin$pred))
lin.df <- ts(lin.df$VALUE, frequency = 4, start = c(2011, 1))
tail(lin.df)

# Ok Linear is done, repeat for quadratic now.

######################### Quadratic Model QuadMod orignial, bigger
predict_quad <- predict(quadmod, n.ahead = 5)
predict_quad$pred

quad1.df <- data.frame(VALUE = c(quadmod$fitted, predict_quad$pred))
quad1.df <- ts(quad1.df$VALUE, frequency = 4, start = c(2011, 1))
tail(quad1.df)
# Ok Quad is done, repeat for other quad now.

######################### Quadratic Model Reduced
predict_quad <- predict(quadmod2, n.ahead = 5)
predict_quad$pred

quad.df <- data.frame(VALUE = c(quadmod2$fitted, predict_quad$pred))
quad.df <- ts(quad.df$VALUE, frequency = 4, start = c(2011, 1))
tail(quad.df)
# Ok Quad is done, repeat for Log now.

############################## Log Model
predict_log <- predict(logmod, n.ahead = 5)
predict_log$pred

log.df <- data.frame(VALUE = c(exp(logmod$fitted), exp(predict_log$pred)))
log.df <- ts(log.df$VALUE, frequency = 4, start = c(2011, 1))
tail(log.df)

######################################################
# Now to plot all of them together with new dfs


plot(incidence_ts, col = "black", lty = 2,
     main = "TB Herd Incidence (%) Time Series Models",
     xlab = "Year", ylab = "Herd Incidence (%)", ylim = c(2,5),
     xlim = c(2020,2024))

# Adding the good Quadratic Model
lines(quad1.df, col = "darkorange2")
lines(quad.df, col = "firebrick3")

# Adding the Linear Model(Linear and Logistic proved similar)
lines(lin.df, col = "royalblue2")
# Adding log model?


legend("bottomright", legend = c("Quad Trend no AR", 
                                 "Linear Trend","Quad Trend with AR" ), 
       col = c("firebrick3",  "royalblue2","darkorange2"), lty = 1, 
       cex = 0.8, y.intersp = 1.5, bty = "n")

#lines(log.df, col = "purple3")

plot(forecast(quadmod, h = 5))
plot(forecast(quadmod2, h = 5), ylim = c(2,5.3))
plot(forecast(linmod, h = 5), ylim = c(2,5.3),
     xlim = c(2020, 2024), main = "Final Model Predictions", xlab = "Year",
     ylab = "Herd Incidence %")
lines(incidence_ts, lty = 2)
lines(forecast(quadmod2, h = 5))
plot(forecast(logmod, h = 5))

