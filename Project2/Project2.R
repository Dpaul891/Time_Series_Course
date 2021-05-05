library("quantmod")
library(TSA)
library(forecast)
library(stats)
library(tseries)
library("depmixS4")

getSymbols("BABA")
BABA.all.ts = window(BABA, start="2015-01-01", end="2019-12-31")
BABA.ts = window(BABA$BABA.Close, start="2015-01-01", end="2019-12-31")
plot(BABA.ts)

#Q2
BABA.arima = auto.arima(BABA.ts)

BABA.all.ts$fitted = BABA.all.ts$BABA.Close
for (i in 1:length(BABA.ts)){
  BABA.all.ts$fitted[i] = BABA.arima$fitted[i]
}

plot(BABA.ts)
lines(BABA.all.ts$fitted, col="red", lwd=1.5, lty=2)
checkresiduals(BABA.arima)

BABA.arima.forecast = forecast(BABA.arima, h=21)
plot(BABA.arima.forecast)
BABA.2020.1 = window(BABA$BABA.Close, start="2020-01-01", end="2020-01-31")
BABA.2020.1.forecast = BABA.2020.1
for (i in 1:length(BABA.2020.1)) {
  BABA.2020.1.forecast[i] = BABA.arima.forecast$mean[i]
}
plot(BABA.2020.1)
lines(BABA.2020.1.forecast, col="red", lwd=1.5, lty=2)


##Q3
BABA.arima.res = BABA.arima$residuals
plot((BABA.arima.res-mean(BABA.arima.res))^2)
acf((BABA.arima.res-mean(BABA.arima.res))^2)

BABA.arima.res.garch = garch(BABA.arima.res, order=c(1,1), grad="numerical", trace=FALSE)
summary(BABA.arima.res.garch)
plot(BABA.arima.res)
lines(fitted(BABA.arima.res.garch)[,1], col="red")
lines(fitted(BABA.arima.res.garch)[,2], col="red")
plot(fitted(BABA.arima.res.garch)[,1]^2)

checkresiduals(resid(BABA.arima.res.garch))
qqnorm(resid(BABA.arima.res.garch))
abline(0,1)

#Q5
BABA.return = diff(BABA.ts)
plot(BABA.return)
returns = as.numeric(BABA.return)
hmm = depmix(returns~1, family=gaussian(), nstates=2, data=data.frame(returns=returns))
hmmfit = fit(hmm, verbose=FALSE)
post_probs = posterior(hmmfit)
plot(2-post_probs[,1], type= 'l', lwd=2, ylab='Viterbi decoding')
plot(post_probs[,2], type= 'l', ylab='posterior decoding')

BABA.ts.log = log(BABA.ts)
plot(BABA.ts.log)
BABA.log.return = diff(BABA.ts.log)
plot(BABA.log.return)
returns = as.numeric(BABA.log.return)
hmm = depmix(returns~1, family=gaussian(), nstates=2, data=data.frame(returns=returns))
hmmfit = fit(hmm, verbose=FALSE)
post_probs = posterior(hmmfit)
plot(2-post_probs[,1], type= 'l', lwd=2, ylab='Viterbi decoding')
plot(post_probs[,2], type= 'l', ylab='posterior decoding')
