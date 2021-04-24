library(TSA)
library(stats)
library(forecast)
library(tseries)

#Question 1
data("airpass")
airpass
log.airpass = log(airpass)
plot(airpass)
plot(log.airpass)

first.diff = diff(log.airpass)
plot(first.diff)

season.first.diff = diff(first.diff, 12)
plot(season.first.diff)

acf(season.first.diff)

log.airpass.arima = stats::arima(log.airpass, order=c(0,1,1), seasonal=list(order=c(0,1,1), period=12))
log.airpass.arima.fit = log.airpass - log.airpass.arima$residuals
plot(log.airpass)
lines(log.airpass.arima.fit, col="red", lty=2)

checkresiduals(log.airpass.arima$residuals)
Box.test(log.airpass.arima$residuals, lag=10, type="Ljung", fitdf=1)
shapiro.test(log.airpass.arima$residuals)
log.airpass.arima.forecast = forecast(log.airpass.arima, h=2*12)
plot(log.airpass.arima.forecast)

#Question 4
set.seed(2)
ma.1 = arima.sim(n=10000, list(ma=1))
plot(ma.1)
acf(ma.1)

set.seed(2)
ma.2 = arima.sim(n=10000, list(ma=c(1, 0.5)))
plot(ma.2)
acf(ma.2)

set.seed(2)
ar.1 = arima.sim(n=10000, list(ar=0.5))
plot(ar.1)
acf(ar.1)

set.seed(2)
arma11 = arima.sim(n=10000, list(ar=0.7, ma=0.5))
plot(arma11)
acf(arma11)

#Question 5
data("google")
plot(google)
acf(google)
plot((google-mean(google))^2)
acf((google-mean(google))^2)

google.garch = garch(google, order=c(1,1), grad="numerical", trace=FALSE)
summary(google.garch)
plot(google)
lines(fitted(google.garch)[,1], col="red")
lines(fitted(google.garch)[,2], col="red")
plot(fitted(google.garch)[,1]^2)

shapiro.test(resid(google.garch))
acf(resid(google.garch)[!is.na(resid(google.garch))])
qqnorm(resid(google.garch))
abline(0,1)
checkresiduals(resid(google.garch))

#Question 6
set.seed(2)
ar.1 = arima.sim(n=1000, list(ar=0.7))
plot(ar.1)
acf(ar.1)
ARMAspec(model=list(ar=0.7))

set.seed(2)
ar.1 = arima.sim(n=1000, list(ar=-0.4))
plot(ar.1)
acf(ar.1)
ARMAspec(model=list(ar=-0.4))

#Question 7
set.seed(2)
ma.1 = arima.sim(n=1000, list(ma=-0.6))
plot(ma.1)
acf(ma.1)
ARMAspec(model=list(ma=-0.6))

set.seed(2)
ma.1 = arima.sim(n=1000, list(ma=0.8))
plot(ma.1)
acf(ma.1)
ARMAspec(model=list(ma=0.8))

