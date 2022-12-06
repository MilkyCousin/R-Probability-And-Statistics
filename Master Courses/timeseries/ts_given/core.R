library("forecast")
library("tseries")

# Зчитуємо дані
data <- read.csv('data.csv', skip=14)
data.subset <- data[,c("date", "open")]
date.split <- strsplit(data.subset$date, split='-')

# Виокремлюємо рік, місяць і день
data.subset["year"] <- unlist(lapply(date.split, function(u) u[1]))
data.subset["month"] <- unlist(lapply(date.split, function(u) u[2]))
data.subset["day"] <- unlist(lapply(date.split, function(u) u[3]))

date.split <- NULL

# Збираємо щомісячні дані (перший відомий день кожного місяця)
years <- unique(data.subset$year)
months <- unique(data.subset$month)

len.months <- length(months) 
len.years <- length(years)
len.vals <- len.years * len.months
values <- numeric(len.vals)
for(i in 1:len.years)
{
  cur.year <- years[i]
  for(j in 1:len.months)
  {
    cur.month <- months[j]
    values[(i - 1) * len.months + j] <- data.subset[
      (data.subset$year == cur.year) & (data.subset$month == cur.month),
      ]$open[1]
  }
}
values <- as.numeric(na.omit(values))

# Остаточний часовий ряд
data.ts <- ts(
  data=values, start=c(1970, 1), end=c(2022, 11), frequency=12
)
plot(data.ts, main="Часовий ряд", ylab="Price")
grid()

# Прологарифмуємо
log.data.ts <- log(data.ts)
plot(log.data.ts, main="Часовий ряд після логарифмування", ylab="Price")
grid()

# Цікавлять дані за 2002-2012 роки
sampled.ts <- window(log.data.ts, start=c(2002, 1), end=c(2012, 12))
plot(sampled.ts, main="Вибірка часового рядку", ylab="Price")
grid()

# Спробумо прибрати тренд
delta.sampled.ts <- diff(sampled.ts)
plot(delta.sampled.ts, main="Вибірка часового рядку", ylab="Price")
grid()
abline(h=0, lty=2, col="red")

# Диференціювання дало можливість позбутися тренду. Виходить щось схоже на шум
# Подивимося на автокореляцію процесу та часткову автокореляцію
par(mfrow=c(1,2))
acf(delta.sampled.ts)
pacf(delta.sampled.ts)
par(mfrow=c(1,1))

# 2 / sqrt(length(sampled.ts))

# Є підозора, що диференційований процес є шумом. Застосуємо тест Ljung-Box
p.values.lb <- numeric(20)
for(h in 1:20)
{
  lb.test <- Box.test(delta.sampled.ts, lag=h, type="Ljung-Box")
  p.values.lb[h] <- lb.test$p.value
}
print(p.values.lb)

# Побудова QQ-діаграми для приростів відносно гауссового розподілу
qqnorm(delta.sampled.ts)
qqline(delta.sampled.ts)
grid()

# Перевірка на гауссовість
sw.test <- shapiro.test(delta.sampled.ts[abs(delta.sampled.ts) <= 0.12])
print(sw.test)

# Чи є стаціонарність відносно тренду? Використаємо KPSS тест
KPSS.test <- kpss.test(sampled.ts, null="Trend")
print(KPSS.test)

# Чи є сенс підігнати процес Xt = log(Yt) як процес AR(1)?
# Тобто Xt = phi * X{t-1} + eps{t}
# Моєливо варто підібрати якусь іншу модель, окрім AR(1)?
# Використаємо критерій Акаїке, робимо мінімізацію по p, q.

# Автокореляція та частинна автокореляція процесу
par(mfrow=c(1,2))
acf(sampled.ts)
pacf(sampled.ts)
par(mfrow=c(1,1))

# Перевіримо умову на коефіцієнт phi за допомогою тесту Dickey-Fuller

# Augmented Dickey-Fuller test
df.test <- adf.test(sampled.ts)
print(df.test)

# Діаграма розсіювання майбутніх значень на попередні
plot(sampled.ts[-1], sampled.ts[-length(sampled.ts)], 
     xlab="X_{t-1}", ylab="X_t",
     main="Лінійна форма залежності майбутнього від минулого")
grid()

# Тобто ситуація близиться до того, що коефіцієнт phi близький до одиниці.
# Отже, процес може представляти з себе випадкове блукання.

# Моделювання AR(1) = ARIMA(1,0,0)
model.ar1 <- arima(sampled.ts, order=c(1,0,0))
print(model.ar1)

# Моделювання MA(1) = ARIMA(0,0,1)
model.ma1 <- arima(sampled.ts, order=c(0,0,1))
print(model.ma1)

# Моделювання I(1) = ARIMA(0,1,0)
model.rw1 <- arima(sampled.ts, order=c(0,1,0))
print(model.rw1)

# Прогнозування

# Зробимо прогноз ціни на рік вперед
years.next <- 1

par(mfrow=c(1,2))

# ARIMA(1,0,0)
Arima.model <- Arima(sampled.ts, order=c(1,0,0))
pred.test <- forecast(Arima.model, h=years.next * 12)
true.test <- window(log.data.ts, start=c(2013, 1), end=c(2013, 12))

plot(pred.test, ylim=c(min(sampled.ts), max(pred.test$mean, true.test)))
lines(window(log.data.ts, start=c(2012, 12), end=c(2013, 12)), lty=2)
grid()

# ARIMA(0,1,0)
Arima.model <- Arima(sampled.ts, order=c(0,1,0))
pred.test <- forecast(Arima.model, h=years.next * 12)
true.test <- window(log.data.ts, start=c(2013, 1), end=c(2013, 12))

plot(pred.test, ylim=c(min(sampled.ts), max(pred.test$mean, true.test)))
lines(window(log.data.ts, start=c(2012, 12), end=c(2013, 12)), lty=2)
grid()

par(mfrow=c(1,1))

# Моделювання на рік вперед

#set.seed(121)
model.test <- c(sampled.ts[length(sampled.ts)], numeric(length(pred.test)))
eps.test <- rnorm(length(pred.test), mean=0, sd=sqrt(model.rw1$sigma2))
for(j in 1:length(pred.test))
{
  model.test[j+1] <- model.test[j] + eps.test[j]
}
model.test <- ts(model.test, start=c(2012, 12), end=c(2013, 12), frequency=12)

#plot(window(log.data.ts, start=c(2010, 1), end=c(2013, 12)), ylim=c(
#  min(sampled.ts, model.test, true.test), 
#  max(sampled.ts, model.test, true.test)
#), ylab="Price", main="Modeled future values using ARIMA(0,1,0)")
lines(model.test, col="darkblue", lty=2)
grid()
