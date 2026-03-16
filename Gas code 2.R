#install.packages("dplyr")
#install.packages("tseries")
install.packages("forecast")
no#install.packages("genetics")

library(tseries)
library(forecast)
library(genetics)
library(dplyr)


gas = read.csv(file = "nrg_cb_gasm__custom_9879805_linear.csv", header = TRUE, sep = ",")
data = gas[c(1:192),]
data$month = as.Date(paste(data$TIME_PERIOD, "-01", sep=""))


df = data %>% 
  filter(geo == "DE") %>% 
  dplyr::select(month, OBS_VALUE)

colnames(df) = c("time","X")
plot(df$time,df$X, type = "l")
par(mfrow=c(1,2))
acf(df$X)
pacf(df$X)
#if there is spikes it may mean it is not a white noice process
n = 192
month = c(rep(seq(1:12),floor(n/12)))

fmonth = as.factor(month)

mod1 = lm(df$X~0+fmonth)

plot(predict(mod1,df$time),type = "l")
plot(df$time,mod1$residuals,type = "l")
acf(mod1$residual)
pacf(mod1$residual)


xc = array(data = 0, dim = c(6,n))
xs = array(data = 0, dim = c(5,n))

f = (1:5)/12
for (t in 1:n) {
  for (j in 1:5) {
    xc[j,t]= cos(2*pi*f[j]*t)
    xs[j,t]= sin(2*pi*f[j]*t)
  }
  xc[6,t] =(-1)^t
}

Xc = as.data.frame(t(xc))
Xs = as.data.frame(t(xs))

par(mfrow=c(1,2))
matplot(Xc[1:3],type = "l")
matplot(Xs[1:3],type = "l")

X = cbind(Xc,Xs)
colnames(X) = c("c1","c2","c3","c4","c5","c6","s1","s2","s3","s4","s5")

mod2 <- lm(df$X~., data = X)
#summary(mod2)

mod3 <- lm(df$X~c1+s1, data = X)
plot(predict(mod3,X),type = "l")
lines(predict(mod2,X),col="blue")


plot(mod3$residuals,type = "l")
acf(mod3$residuals)


mod4 = arima(df$X,order = c(12,0,0), include.mean=FALSE)
par(mfrow=c(1,3))
plot(mod4$residuals, type = "l")
acf(mod4$residuals, main="")
pacf(mod4$residuals, main="")

res = mod4$residuals
res.acf = acf(res, plot=FALSE)
Lb.result = Box.test(res, lag = length(res.acf$acf)-1,
                     type = "Ljung-Box", fitdf = 12)
#p value is over 0.05 so the residuals may be from a white noise process

coef(mod4)
mod4

mod5 = arima(df$X, order = c(12,1,0), fixed=c(NA,NA,0,0,0,0,0,0,0,0,NA,0), include.mean=FALSE)
mod5
#mod4 smaller aic

xlog = log(df$X)
par(mfrow=c(1,2))
plot(xlog, type = "l")
plot(df$X, type = "l")

acf(xlog, main="")
pacf(xlog, main="")

dx = diff(xlog)
acf(dx)
pacf(dx)

mod6 = arima(dx, order = c(9,0,0), include.mean = FALSE)
par(mfrow=c(3,1))
r = mod6$residuals
plot(r)
acf(r,main="")
pacf(r,main="")

r.acf = acf(r, plot=FALSE)
Lb.result1 = Box.test(r, lag = length(res.acf$acf)-1,
                     type = "Ljung-Box", fitdf = 9)
#it looks reasonable
mod6$coef

pacf.theory = ARMAacf(ar = 0, ma=mod6$coef, lag.max = 30, pacf=TRUE)
acf.theory = ARMAacf(ar = 0, ma=mod6$coef, lag.max = 30, pacf=FALSE)

par(mfrow=c(1,2))
acf(dx, main="")
lines(0:30,acf.theory, col="red")
pacf(dx, main="")
lines(pacf.theory, col="red")

#forecasting
m = 19
n1 = length(xlog)-m
mindex = (n1+1):(n1+m)
xtrain = xlog[1:n1]
xtest = xlog[mindex]

par(mfrow=c(1,1))
plot(1:n,xtrain[1:n],xlim = c(1,n+m),col="blue",type = "l")
points(mindex,xtest,col="red",type = "l")

xf = 