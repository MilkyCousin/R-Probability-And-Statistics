## Підготовка даних до роботи
data <- read.csv('city_temperature.csv')
print(sort(unique(data$Country)))

data.ukraine <- data[data$Country == "Ukraine", ]
print(unique(data.ukraine$City))

data.kiev <- data.ukraine[
  data.ukraine$City == "Kiev", c("Day", "Month", "Year", "AvgTemperature")
  ]

min.year <- min(data.kiev$Year)
min.month <- min(data.kiev$Month)

max.year <- max(data.kiev$Year)
max.month <- max(data.kiev$Month)

plot(data.kiev$AvgTemperature, type='l')

#print(length(data.kiev$AvgTemperature))
data.kiev.ts <- ts(
  data=data.kiev$AvgTemperature, 
  start=c(min.year, min.month), 
  end=c(max.year, max.month),
  frequency=365
  )
#print(length(data.kiev.ts))
data.kiev.ts <- (data.kiev.ts - 32) * 5/9

plot(data.kiev.ts,
     ylab="Temperature",
     main="Щоденна середня температура у м.Київ")
grid()

sampled.ts <- window(data.kiev.ts, 
                     start=c(2005, 1), 
                     end=c(2010, 12),
                     frequency=365)
n <- length(sampled.ts)

plot(sampled.ts,
     ylab="Temperature",
     main="Щоденна середня температура у м.Київ, 2005-2010")
grid()

out.mask <- abs(sampled.ts - mean(sampled.ts)) > 4 * sd(sampled.ts)
print(sampled.ts[out.mask])

abline(h=mean(sampled.ts) - 4 * sd(sampled.ts), col='red', lty=2)

# Заміна локальним середнім
idx <- (1:n)
eps <- 7
for(j in idx[out.mask])
{
  x0 <- sampled.ts[j]
  eps.neighborhood <- abs(idx[-j] - j) < eps
  sampled.ts[j] <- sum(eps.neighborhood * sampled.ts[-j]) / sum(eps.neighborhood[-j])
}
print(sampled.ts[out.mask])

plot(sampled.ts,
     ylab="Temperature",
     main="Щоденна середня температура у м.Київ, 2005-2010")
grid()

## Аналіз даних для вибору моделі

### Вгадування тренду, підгонка параметрів за допомогою МНК

poly.mean <- function(t, A, B, C, D)
{
  t0 <- t - 2005
  A + B * t0 + C * t0^2 + D * t0^3
}

sampled.values <- as.numeric(sampled.ts)
sampled.time <- time(sampled.ts)
data.ts <- data.frame(cbind(sampled.time, sampled.ts))
colnames(data.ts) <- c("X", "Y")

poly.ols <- lm(Y ~ I(X - 2005) + I((X - 2005)^2) + I((X - 2005)^3),
               data=data.ts)
coef.poly.ols <- coef(poly.ols)

poly.mean.ols <- function(t) 
{ 
  poly.mean(t, 
            A=coef.poly.ols[1], 
            B=coef.poly.ols[2], 
            C=coef.poly.ols[3], 
            D=coef.poly.ols[4]) 
}
curve(poly.mean.ols(x), col="red", lty=2, add=T)
sampled.ts.const <- sampled.ts - poly.mean.ols(sampled.time)
data.ts$Y1 <- as.numeric(sampled.ts.const)

plot(sampled.ts.const,
     ylab="Temperature",
     main="Polynomial centered time series")
grid()

trig.mean <- function(t, A, B, C, D)
{
  A * sin(B * t + C) + D
}

trig.ols <- nls(Y1 ~ A * sin(B * X + C) + D, 
                data=data.ts,
                start=list(A=15, B=2*pi, C=4.5, D=-3))
coef.trig.ols <- coef(trig.ols)

trig.mean.ols <- function(t)
{
  trig.mean(t, 
            A=coef.trig.ols["A"], 
            B=coef.trig.ols["B"], 
            C=coef.trig.ols["C"],
            D=coef.trig.ols["D"])
}

curve(trig.mean.ols(x), col="red", lty=2, lwd=2, add=T)
# what the fuck?

data.ts$Y2 <- data.ts$Y1 - trig.mean.ols(data.ts$X)
plot(data.ts$X, data.ts$Y2, type="l")

full.nls <- nls(Y ~ A * sin(B * X) + C + D * (X - 2005) + E * (X - 2005)^2,
                data=data.ts,
                start=list(
                  A=15, B=2*pi, C=5.32, 
                  D=4.37, E=-1.28
                  )
                )
summary(full.nls)

full.func <- function(t)
{
  A <- coef(full.nls)['A']
  B <- coef(full.nls)['B']
  C <- coef(full.nls)['C']
  D <- coef(full.nls)['D']
  E <- coef(full.nls)['E']
  A * sin(B * t) + C + D * (t - 2005) + E * (t - 2005)^2
}

plot(sampled.ts)
curve(full.func(x), col="red", lty=2, lwd=2, add=T)
grid()

res.nls <- sampled.ts - fitted(full.nls)

par(mfrow=c(1,3))

plot(res.nls, main="Residuals")
grid()

qqnorm(res.nls)
qqline(res.nls)
grid()

hist(res.nls, prob=T)
curve(dnorm(x, mean(res.nls), sd(res.nls)), col='red', add=T)
grid()

par(mfrow=c(1,1))

plot(fitted(full.nls), resid(full.nls), cex=0.25, xlab="Прогноз", ylab="Залишки")

# Перевірка на стаціонарність процесу залишків

par(mfrow=c(1,2))
acf(res.nls)
pacf(res.nls)
par(mfrow=c(1,1))

p.max <- 14
q.max <- 3

AIC.matrix <- matrix(nrow=p.max+1, ncol=q.max+1)

for(p in 0:p.max)
{
  for(q in 0:q.max)
  {
    if((0 <= p + q) & (p + q <= min(p.max, q.max)))
    {
      arma.pq <- arima(res.nls, order=c(p, 0, q))
      AIC.matrix[p+1,q+1] <- AIC(arma.pq)
    }
  }
}

AIC.matrix

res.nls.arima <- arima(res.nls, order=c(3,0,0))

library("tseries")
kpss.test.res <- kpss.test(res.nls)
kpss.test.res

## Прогноз погоди

to.forecast <- window(data.kiev.ts, 
                      start=c(2011, 1), 
                      end=c(2012, 1),
                      frequency=365)
times.for.forecast <- time(to.forecast)
data.ts.new <- data.frame(cbind(times.for.forecast, as.numeric(to.forecast)))
colnames(data.ts.new) <- c("X", "Y")

par(mfrow=c(1,1))
predicted.temp <- predict(full.nls, data.ts.new)
plot(times.for.forecast, as.numeric(to.forecast), 
     type='l', xlab='Days', ylab='Temperature', main='Forecast for 2011')
curve(full.func(x), col='red', add=T)
grid()

plot(full.func(times.for.forecast) - as.numeric(to.forecast),  type='l',
     xlab='Days', ylab='Residuals', main='Forecast for 2011')
abline(h=0, lty=2, col='red')
grid()
