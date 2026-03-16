library(dplyr)
library(tseries)

#load data
data = read.csv("nrg_cb_em_linear.csv")
data = filter(data, siec=="E7000" & unit=="GWH" & nrg_bal=="AIM")


data$month = as.Date(paste(data$TIME_PERIOD, "-01", sep=""))

df = data %>% filter(geo =="DE") %>% select(month,OBS_VALUE)
colnames(df) = c("time", "X")

min(df$time)
max(df$time)

plot(df$time, df$X, type="l")#raw original data

print(df)
summary(df$X)
boxplot(df$X)


xlog = log(df$X)#log of data
plot(xlog, type="l")

acf(xlog,main="")#check if white noise
pacf(xlog,main="")

#check stationarity

adf.test(xlog, alternative = "stationary", k = 13)

#seasonal difference
sxlog = diff(xlog, lag=12)

adf.test(sxlog, alternative = "stationary", k = 13)

#first and seasonal difference
dsxlog = diff(sxlog)

adf.test(dsxlog, alternative = "stationary", k = 13)#is now stationary

plot(dsxlog,type="l")
acf(dsxlog,main="")
pacf(dsxlog,main="")

#suggests seasonal and first order differencing is appropriate

#The significant spike at lag 1 and 2 in the ACF suggests a non-seasonal MA(1) component or a non-seasonal MA(2) component
#the significant spike at lag 12 in the ACF suggests a seasonal MA(1) component 
#begin with an ARIMA(0,1,2)(0,1,1) model, indicating a first and seasonal difference, and non-seasonal MA(2) and seasonal MA(1) components.

mod1 = arima(xlog, order = c(0, 1, 2), seasonal=list(order=c(0,1,1), period=12))

r = mod1$residuals
BT = Box.test(r, lag=(length(acf(r))-1), type = "Ljung-Box", fitdf=1)

acf(mod1$residuals,main="")
pacf(mod1$residuals,main="")
AIC(mod1)

#significant at lag 6 and almost significant at lag 14
#suggests more non-seasonal terms need to be added

#trying different non-seasonal changes
mod00 = arima(xlog, order = c(0, 1, 0), seasonal=list(order=c(0,1,1), period=12))
mod00
mod01 = arima(xlog, order = c(0, 1, 1), seasonal=list(order=c(0,1,1), period=12))
mod01
mod02 = arima(xlog, order = c(0, 1, 2), seasonal=list(order=c(0,1,1), period=12))
mod02
mod03 = arima(xlog, order = c(0, 1, 3), seasonal=list(order=c(0,1,1), period=12))
mod03
mod10 = arima(xlog, order = c(1, 1, 0), seasonal=list(order=c(0,1,1), period=12))
mod10
mod20 = arima(xlog, order = c(2, 1, 0), seasonal=list(order=c(0,1,1), period=12))
mod20
mod30 = arima(xlog, order = c(3, 1, 0), seasonal=list(order=c(0,1,1), period=12))
mod30
mod11 = arima(xlog, order = c(1, 1, 1), seasonal=list(order=c(0,1,1), period=12))
mod11
mod22 = arima(xlog, order = c(2, 1, 2), seasonal=list(order=c(0,1,1), period=12))
mod22
mod33 = arima(xlog, order = c(3, 1, 3), seasonal=list(order=c(0,1,1), period=12))
mod33
mod12 = arima(xlog, order = c(1, 1, 2), seasonal=list(order=c(0,1,1), period=12))
mod12
mod13 = arima(xlog, order = c(1, 1, 3), seasonal=list(order=c(0,1,1), period=12))
mod13
mod21 = arima(xlog, order = c(2, 1, 1), seasonal=list(order=c(0,1,1), period=12))
mod21
mod23 = arima(xlog, order = c(2, 1, 3), seasonal=list(order=c(0,1,1), period=12))
mod23
mod31 = arima(xlog, order = c(3, 1, 1), seasonal=list(order=c(0,1,1), period=12))
mod31
mod32 = arima(xlog, order = c(3, 1, 2), seasonal=list(order=c(0,1,1), period=12))
mod32

AIC(mod10)

#aic of arima(1,1,0)(0,1,1) = -732.3188 is lowest

tsdiag(mod10)
r = mod10$residuals
BT = Box.test(r, lag=(length(acf(r))-1), type = "Ljung-Box", fitdf=1)
BT

#Looking at the ACF/PACF of the residuals appears to be consistent with white noise. LB test is conclusive. 
#The ARIMA(1,1,0)(0,1,1) model for xlog appears reasonable.

#forecasting
library(forecast)
#10%=19 20%=38
m = 19
n = length(xlog) - m
mindex = (n+1):(n+m)
Xtrain = xlog[1:n]
Xtest = xlog[mindex]

plot(1:n, Xtrain[1:n], xlim=c(1,n+m), col="blue",type="l")#shows split data
points(mindex, Xtest,col="red", type="l")

mod = arima(Xtrain, order = c(1, 1, 0), seasonal=list(order=c(0,1,1), period=12))
mod.f = predict(mod, n.ahead = m)


xf1 = mod.f$pred

#95% confidence interval

xf.l = mod.f$pred - 1.96*mod.f$se
xf.u = mod.f$pred + 1.96*mod.f$se

#plot on same graph

ylims = c(min(df$X),max(df$X))
plot(1:n, Xtrain[1:n], xlim = c(1,n+m), col="blue",type = "l")
points(mindex, Xtest,col="red", type="l")
lines(mindex, xf1,col="green")
lines(mindex, xf.l,col="green", lty=2)
lines(mindex, xf.u,col="green", lty=2)


#adjust for a constant after differencing

mod2 = forecast::Arima(ts(Xtrain), order = c(1, 1, 0), seasonal=list(order=c(0,1,1), period=12), include.constant = TRUE)
mod.f2 = forecast::forecast(mod2, h=m, level = c(90, 95))
xf2 = mod.f2$mean
# We use 95% interval
xf2.l = mod.f2$lower[,2]
xf2.u = mod.f2$upper[,2]

#plot
ylims = c(min(df$X),max(df$X))
plot(1:n, Xtrain[1:n], xlim = c(1,n+m), col="blue",type = "l")
points(mindex, Xtest,col="red", type="l")
lines(mindex, xf1,col="green")
lines(mindex, xf2.l,col="green", lty=2)
lines(mindex, xf2.u,col="green", lty=2)


#forcast evaluation
#rolling origin, no model update

#store differences and forecast
d = vector()
xf = vector()

#evaluate forecast on test set
for (i in 1:m){
  # update info
  x.this = xlog[1:n+i-1]
  mod.this = arima(x.this, order = c(1, 1, 0), seasonal=list(order=c(0,1,1), period=12), fixed = as.numeric(mod$coef))
  mod.f = predict(mod.this, n.ahead=1)
  xf[i] = mod.f$pred
  d[i] = xf[i]-xlog[n+i]
  
}

#plot
ylims = c(min(df$X),max(df$X))
plot(1:n, Xtrain[1:n], xlim = c(1,n+m), col="blue",type = "l")
points(mindex, Xtest,col="red", type="l")
lines(mindex[1:m], xf,col="green")
lines(mindex, xf2.l,col="green", lty=2)
lines(mindex, xf2.u,col="green", lty=2)

d = vector()
xf = vector()
k = 2
#evaluate forecast on test set
for (i in 1:(m-(k-1))){
  # update info
  x.this = xlog[1:n+i-1]
  mod.this = arima(x.this, order = c(1, 1, 0), seasonal=list(order=c(0,1,1), period=12), fixed = as.numeric(mod$coef))
  mod.f = predict(mod.this, n.ahead=k)
  xf[i] = mod.f$pred[k]
  d[i] = xf[i]-xlog[n+i-1+k]
  
}

#plot
ylims = c(min(df$X),max(df$X))
plot(1:n, Xtrain[1:n], xlim = c(1,n+m), col="blue",type = "l")
points(mindex, Xtest,col="red", type="l")
lines(mindex[k:m], xf,col="green")
lines(mindex, xf2.l,col="green", lty=2)
lines(mindex, xf2.u,col="green", lty=2)

#Calculate error metrics

rmse = function(d){
  return(mean(abs(d^2)))
}

mae = function(d){
  return(mean(abs(d)))  
}

mape = function(d,xlog){
  return(mean(100*abs(d/xlog)))  
}

mape(d,Xtest[k:m])

futp = predict(mod10,6)
L95 = exp(futp$pred - 1.96*futp$se)
U95 = exp(futp$pred + 1.96*futp$se)

Forecast = exp(futp$pred + mod$sigma2/2)

dff = data.frame(L95,Forecast,U95)
print(dff,row.names=FALSE)