# Install.packages("dplyr")
# Install.packages("tseries")
# Install.packages("forecast")
# Install.packages("generics")
library(dplyr)
library(tseries)
library(forecast)
library(generics)

# Load data into R
Germany_oil <- read.csv("~/Downloads/Uni/Math 334/group project/code/oil data.csv")

# Information we use for model creation formatting
data <- Germany_oil[c(1:191),]
data$month <- as.Date(paste(data$TIME_PERIOD, "-01", sep = ""))
df = data %>% select(month, OBS_VALUE)
colnames(df) <- c("time", "THS_T")
par(mfrow = c(1,3))
plot(df$time, df$THS_T, type = "l")

# checking for white noises process using acf function
acf(df$THS_T) # because there are spikes in lags which are outside of the blue line,
# this indicates that this might not be a white noise process.
pacf(df$THS_T)

# modelling a seasonal trend
n <- length(df$THS_T)
month <- c(rep(seq(1:12), floor(n/12)), 1:(n %% 12))
fmonth <- as.factor(month)
mod1 <- lm(df$THS_T ~ 0 +fmonth)
summary(mod1)
plot(predict(mod1, df$time), type = "l")

# model mod residuals after removing seasonal trend
plot(df$time, mod1$residuals, type = "l")
par(mfrow = c(1,2))
acf(mod1$residuals)
pacf(mod1$residuals)

# The ACF and PACF we get after removing seasonal trend looks weirder than not removing seasonal trend, therefore, we stick to our original data

diffdata <- diff(df$THS_T)
plot(diffdata, type = "l")

modestimate <- arima(diffdata,c(0,0,12), include.mean = FALSE)
acf(modestimate$residuals)
pacf(modestimate$residuals)
tsdiag((modestimate))

summary(modestimate)

#need to use coeff function to check which coefficient needs to be restricted.
modestimate1 <- arima(diffdata, c(0,0,12), fixed = c(NA,NA,NA,NA,0,NA,0,NA,0,NA,0,NA), include.mean = FALSE)
acf(modestimate1$residuals)
pacf(modestimate1$residuals)
tsdiag(modestimate1)

main_model <-  arima(diffdata[1:179], c(0,0,12), fixed = c(NA,NA,NA,NA,0,NA,0,NA,0,NA,0,NA), include.mean = FALSE)
forecasting_testdata <- predict(main_model, n.ahead = 12)

xf = forecasting_testdata$pred
xf.l  = xf - 1.96*forecasting_testdata$se
xf.u = xf + 1.96*forecasting_testdata$se

par(mfrow= c(1,1))
plot(as.numeric(diffdata[180:191]), type = "l", ylim = c(-1000,1000))
lines(as.numeric(xf), col = "red")
lines(as.numeric(xf.l), col = "green")
lines(as.numeric(xf.u), col = "green")

#Forecasting
main_model2 <- forecast::Arima(ts(diffdata[1:179]), order=c(0,0,12), include.constant = TRUE)
forecasting_testdata2 <- forecast::forecast(main_model2, h=12, level = c(90,95))
xf2 <- forecasting_testdata2$mean
#we use 95% interval
xf2.l <- forecasting_testdata2$lower[,2]
xf2.u <- forecasting_testdata2$upper[,2]
#one way to plot
par(mfrow = c(1,1))
plot(forecasting_testdata2, main= "", ylab = "log THS_T", xlab = "time")

#store differences and forecast
d = vector()
xf = vector()

#evaluate forecast on test set no rolling origin
for (i in 1:12){
  # update info
  x.this = diffdata[1:179+i-1]
  mod.this = arima(x.this, order = c(0,0,12), seasonal=list(order=c(0,0,12), period=12, fixed = as.numeric(forecasting_testdata$coef)))
  mod.f = predict(mod.this, n.ahead=1)
  xf[i] = forecasting_testdata$pred
  d[i] = xf[i]-diffdata[179+i]
}

# evaluating forecast on test set with rolling origin
d = vector()
xf = vector()
k = 12
#evaluate forecast on test set
for (i in 1:(12-(k-1))){
  # update info
  x.this = diffdata[1:179+i-1]
  mod.this = arima(x.this, order = c(0,0,12), seasonal=list(order=c(0,0,12), period=12), fixed = as.numeric(forecasting_testdata$coef))
  mod.f = predict(mod.this, n.ahead=k)
  xf[i] = mod.f$pred[k]
  d[i] = xf[i]-diffdata[179+i-1+k]
}


# Predictions for the upcoming 3 months
Q1model<- arima(diffdata, c(0,0,12), fixed =c(NA,NA,NA,NA,0,NA,0,NA,0,NA,0,NA), include.mean = FALSE)
forecast_Q1 <- predict(Q1model, n.ahead = 3)

Q1 <- forecast_Q1$pred
Q1.l <- Q1 - 1.96* forecast_Q1$se
Q1.U <- Q1 +1.96*forecast_Q1$se

plot(155:191,diffdata[155:191], type = "l", xlim =c(155,194))
lines(192:194, Q1, col ="red")
lines(191:192, c(diffdata[190], Q1[1]), col = "red")

#prediction of main data 
january = df$THS_T[191] + Q1[1]
february = january  + Q1[2]      
march = february + Q1[3]

#plot on the original data
plot(df$THS_T[155:191], type = "l", xlim = c(0,45))
lines(38:40,c(january,february,march), col = "red")
lines(37:38, c(df$THS_T[191],january), col = "red")

# Marrix
rmse = function(d){
  return(mean(abs(d^2)))
}

mae = function(d){
  return(mean(abs(d)))  
}

mape = function(d,xlog){
  return(mean(100*abs(d/diffdata)))  
}

mape(d,155:191)

futp = predict(main_model,6)
L95 = exp(futp$pred - 1.96*futp$se)
U95 = exp(futp$pred + 1.96*futp$se)

Forecast = exp(futp$pred + main_model2$sigma2/2)

dff = data.frame(L95,Forecast,U95)
print(dff,row.names=FALSE)