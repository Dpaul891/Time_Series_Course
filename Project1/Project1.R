library(dlnm)
library(lubridate)
library(forecast)
library(dplyr)
library(TSA)
library(stats)
library(tseries)
library(urca)
library(vars)
Data = chicagoNMMAPS
Data.year.month = group_by(Data,year,month)
Data.year.month.average = summarise(Data.year.month, cvd=mean(cvd, na.rm=T), death=mean(death, na.rm = T), temp=mean(temp, na.rm=T), pm10=mean(pm10,na.rm = T), o3=mean(o3,na.rm = T))

#Question one: Focus on the univariate analysis of the cardiovascular death.
cvd.ts = ts(Data.year.month.average$cvd, start=c(1987,1), freq=12)
cvd.ts.training = window(cvd.ts, start=c(1987,1), end=c(1995,12))
cvd.ts.validation = window(cvd.ts, start=c(1996,1), end=c(1997,12))
cvd.ts.testing = window(cvd.ts, start=c(1998,1), end=c(2000,12))

# mode1 stl
cvd.stl = stl(cvd.ts.training, s.window=4, t.window=12, robust=TRUE)
plot(cvd.stl)
cvd.stl.decom = cvd.stl$time.series
cvd.stl.decom.seasonal = cvd.stl.decom[, 1]
cvd.stl.decom.trend = cvd.stl.decom[, 2]
cvd.stl.decom.remainder = cvd.stl.decom[, 3]
acf(cvd.stl.decom.remainder)
plot(cvd.ts.training, lwd=1)
lines(cvd.stl.decom.trend, col="green", lwd=1.5)
lines(cvd.stl.decom.trend+cvd.stl.decom.seasonal, col="red", lty=2)
cvd.stl.validation.forecast = forecast(cvd.stl, h=2*12)
plot(cvd.stl.validation.forecast)
lines(cvd.ts.validation, lty=2, lwd=1.5)
cvd.stl.validation.accuracy = sqrt(sum((cvd.stl.validation.forecast$mean - cvd.ts.validation)**2)/length(cvd.ts.validation))

# mode2 holt-winter
cvd.hw = HoltWinters(cvd.ts.training, seasonal="add")
plot(cvd.hw$fitted)
plot(cvd.hw)
acf(resid(cvd.hw))
cvd.hw.validation.forecast = forecast(cvd.hw, h=2*12)
plot(cvd.hw.validation.forecast)
lines(cvd.ts.validation, lty=2, lwd=1.5)
cvd.hw.validation.accuracy = sqrt(sum((cvd.hw.validation.forecast$mean - cvd.ts.validation)**2)/length(cvd.ts.validation))

# model3 arima
cvd.arima = auto.arima(cvd.ts.training)
plot(cvd.ts.training)
lines(cvd.arima$fitted, col="red", lwd=1.5, lty=2)
plot(cvd.arima$residuals)
acf(cvd.arima$residuals)
checkresiduals(cvd.arima)
cvd.arima.validation.forecast = forecast(cvd.arima, h=2*12)
plot(cvd.arima.validation.forecast)
lines(cvd.ts.validation, lty=2, lwd=1.5)
cvd.arima.validation.accuracy = sqrt(sum((cvd.ts.validation-cvd.arima.validation.forecast$mean)**2)/length(cvd.ts.validation))

#select the best model
#holt-winter is the best
cvd.ts.train.validation = window(cvd.ts, start=c(1987,1), end=c(1997,12))
cvd.best.hw = HoltWinters(cvd.ts.train.validation, seasonal="add")
plot(cvd.best.hw$fitted)
plot(cvd.best.hw)
acf(resid(cvd.best.hw))
Box.test(resid(cvd.best.hw), type="Ljung")
cvd.best.hw.test.forecast = forecast(cvd.best.hw, h=3*12)
plot(cvd.best.hw.test.forecast)
lines(cvd.ts.testing, lty=2, lwd=1.5)
cvd.best.hw.test.accuracy = sqrt(sum((cvd.best.hw.test.forecast$mean - cvd.ts.testing)**2)/length(cvd.ts.testing))

#Question two: Focus on the multivariate analysis of the cardiovascular death.

#model1 SARIMA with external variable
temp.ts = ts(Data.year.month.average$temp, start=c(1987,1), freq=12)
temp.ts.training = window(temp.ts, start=c(1987,1), end=c(1995,12))
temp.ts.validation = window(temp.ts, start=c(1996,1), end=c(1997,12))
temp.ts.testing = window(temp.ts, start=c(1998,1), end=c(2000,12))

pm10.ts = ts(Data.year.month.average$pm10, start=c(1987,1), freq=12)
pm10.ts.training = window(pm10.ts, start=c(1987,1), end=c(1995,12))
pm10.ts.validation = window(pm10.ts, start=c(1996,1), end=c(1997,12))
pm10.ts.testing = window(pm10.ts, start=c(1998,1), end=c(2000,12))

o3.ts = ts(Data.year.month.average$o3, start=c(1987,1), freq=12)
o3.ts.training = window(o3.ts, start=c(1987,1), end=c(1995,12))
o3.ts.validation = window(o3.ts, start=c(1996,1), end=c(1997,12))
o3.ts.testing = window(o3.ts, start=c(1998,1), end=c(2000,12))

layout(1:4)
plot(cvd.ts.training)
plot(temp.ts.training)
plot(pm10.ts.training)
plot(o3.ts.training)
layout(1:1)

multi.arima.external = auto.arima(cvd.ts.training, xreg=cbind(temp.ts.training, pm10.ts.training, o3.ts.training))
multi.arima.external

plot(cvd.ts.training)
lines(multi.arima.external$fitted, col="red", lwd=1.5, lty=2)
checkresiduals(multi.arima.external)
multi.arima.external.validation.forecast = forecast(multi.arima.external, xreg=cbind(temp.ts.validation, pm10.ts.validation, o3.ts.validation))
plot(multi.arima.external.validation.forecast)
lines(cvd.ts.validation, lty=2, lwd=1.5)
multi.arima.external.validation.accuracy = sqrt(sum((multi.arima.external.validation.forecast$mean - cvd.ts.validation)**2)/length(cvd.ts.validation))
multi.arima.external.validation.accuracy

#model2 vector AR
multi.training = Data.year.month.average[1:length(cvd.ts.training),c("cvd","temp","pm10","o3")]
VARselect(multi.training, lag.max=2, type="const")
cvd.var = VAR(multi.training, p=2, type="const")
plot(multi.training[["cvd"]], type="l")
lines(c(rep(NA, 2), fitted(cvd.var)[,1]), col=2)
checkresiduals(resid(cvd.var)[,"cvd"])
cvd.var
cvd.var.val.farecast = rnorm(length(cvd.ts.validation))
cvd.var.coef = c(0.36411984,-0.29789871,-0.04243912,0.13608390,0.20698752,0.38811791,-0.01432600,-0.28986066,26.80201828)

cvd.ts.train.val = window(cvd.ts, start=c(1987,1), end=c(1997,12))
temp.ts.train.val = window(temp.ts, start=c(1987,1), end=c(1997,12))
pm10.ts.train.val = window(pm10.ts, start=c(1987,1), end=c(1997,12))
o3.ts.train.val = window(o3.ts, start=c(1987,1), end=c(1997,12))

for (ym in (length(cvd.ts.training)+1):(length(cvd.ts.training)+length(cvd.ts.validation))){
  cvd.var.val.farecast[ym-length(cvd.ts.training)] = cvd.var.coef[1]*cvd.ts.train.val[ym-1]+
                                                     cvd.var.coef[2]*temp.ts.train.val[ym-1]+
                                                     cvd.var.coef[3]*pm10.ts.train.val[ym-1]+
                                                     cvd.var.coef[4]*o3.ts.train.val[ym-1]+
                                                     cvd.var.coef[5]*cvd.ts.train.val[ym-2]+
                                                     cvd.var.coef[6]*temp.ts.train.val[ym-2]+
                                                     cvd.var.coef[7]*pm10.ts.train.val[ym-2]+
                                                     cvd.var.coef[8]*o3.ts.train.val[ym-2]+
                                                     cvd.var.coef[9]
  
  cvd.ts.train.val[ym] = cvd.var.val.farecast[ym-length(cvd.ts.training)]
}
cvd.var.val.farecast.ts = ts(cvd.var.val.farecast, start=c(1996,1), freq=12)
plot(cvd.ts.train.validation)
lines(cvd.var.val.farecast.ts, type="l", lwd=1.5, col="blue")
cvd.var.validation.accuracy = sqrt(sum((cvd.var.val.farecast - cvd.ts.validation)**2)/length(cvd.ts.validation))
cvd.var.validation.accuracy
#run the best: ARIMA with external variables

arima.external.best = auto.arima(cvd.ts.train.validation, xreg=cbind(temp.ts.train.val, pm10.ts.train.val, o3.ts.train.val))
arima.external.best
plot(cvd.ts.train.validation)
lines(arima.external.best$fitted, col="red", lwd=1.5, lty=2)

checkresiduals(arima.external.best)

arima.external.best.test.forecast = forecast(arima.external.best, xreg=cbind(temp.ts.testing, pm10.ts.testing, o3.ts.testing))
plot(arima.external.best.test.forecast)
lines(cvd.ts.testing, lty=2, lwd=1.5)
arima.external.best.test.accuracy = sqrt(sum((arima.external.best.test.forecast$mean - cvd.ts.testing)**2)/length(cvd.ts.testing))
arima.external.best.test.accuracy

