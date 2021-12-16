library(rgl)
options(rgl.printRglwidget = TRUE)

path <- "/home/coolpenguin/R/compsta/lab2"

source(paste(path, "qqplotinterval.R", sep="/"))

cur.data <- read.csv(paste(path, "merged.csv", sep="/"),
                     row.names = "X")

lim.val <- 0.000001
rem <- function(val) { ifelse(abs(val) < lim.val, 0, val) }

N.test <- 20
N.last <- 50
# Прогнозуємо ціни закриття компанії cl (Colgate-Palmolive Company)

x.dat.full <- cur.data[-nrow(cur.data),-1]
x.dat.full$cl <- cur.data$cl[-1]
n.full <- nrow(x.dat.full)

# Повні дані, без останніх двадцяти сесій
x.dat <- x.dat.full[1:(n.full-N.test),]
n.dat <- nrow(x.dat)

# 50 сесій перед двадцятьма останніми
x.red <- x.dat.full[(n.full-N.test-N.last+1):(n.full-N.test),]

# Останні двадцять сесій
x.dat.to.predict <- x.dat.full[(n.full-N.test+1):n.full,]

# Аналіз головних компонент за повними даними
pca.full <- princomp(x.dat[,-1], cor = T)
print(summary(pca.full))

# Діаграма власних чисел - злам після першої компоненти, однак
# вагома інформація все ще є для наступних трьох
plot(pca.full)

# Тому спробуємо робити підгонку за першими чотирьма компонентами
P.full <- 4

# Переходимо до першого базису
x.dat.pc <- data.frame(pca.full$scores[,1:P.full])
# Навантаження у нас такі
loadings.full <- apply(pca.full$loadings[], 2, rem)
print(loadings.full[,1:P.full])

lm.full <- lm(x.dat$cl ~ ., data=x.dat.pc)
print(summary(lm.full))
#plot(lm.full)

# Аналіз головних компонент за 50 сесіями
pca.red <- princomp(x.red[,-1], cor = T)
print(summary(pca.red))

# Діаграма власних чисел - злам після першої компоненти, однак
# вагома інформація все ще є для наступних чотирьох
plot(pca.red)

# Тому спробуємо робити підгонку за першими трьома компонентами.
# Було б більше даних і була б така сама картина - зробили б для 5 компонент
P.red <- 4

# Переходимо до другого базису
x.red.pc <- data.frame(pca.red$scores[,1:P.red])

# Навантаження у цьому разі такі
loadings.red <- apply(pca.red$loadings[][], 2, rem)
print(loadings.red[,1:P.red])

lm.red <- lm(x.red$cl ~ ., data=x.red.pc)
print("pcj, j=1,2,3,4")
print(summary(lm.red))

# Значущими є лише зсув та коефіцієнт при 1-ій та 3-ій компонентах. Врахуємо це
lm.red.pc1 <- lm(x.red$cl ~ Comp.1, data=x.red.pc)
print("pc1")
print(summary(lm.red.pc1))

lm.red.pc13 <- lm(x.red$cl ~ Comp.1 + Comp.3, data=x.red.pc)
print("pc1, pc3")
print(summary(lm.red.pc13))

# Чудо-функция. Костыльная

get.original.coef <- function(x, pca.m, pca.lm)
{
  # x - регресори, що використовувалися для princomp(x, cor=T)
  # pca.m <- princomp(x, cor=T)
  # pca.lm - регресія на головні компоненти
  means <- apply(x, 2, mean)
  sdevs <- apply(x, 2, sd)
  m.s <- means/sdevs
  pca.loadings <- pca.m$loadings[]
  pca.coef <- coef(pca.lm)
  names.coef <- names(pca.coef)
  b <- c()
  i <- names.coef == "(Intercept)"
  if(length(i) > 0)
  {
    s.bias <- pca.coef[i]
    s <- 0
    for(name.comp in names.coef[-i])
    {
      print(name.comp)
      p <- colnames(pca.loadings) == name.comp
      clp <- pca.loadings[,p]
      clp.ms <- crossprod(clp, m.s)
      s <- s + clp.ms * pca.coef[names.coef == name.comp]
    }
    a0 <- s.bias - s
    b <- c(b, a0)
    names(b)[1] <- "(Intercept)"
  }
  
  for(name.base in colnames(x))
  {
    s.j <- 0
    l <- rownames(pca.loadings) == name.base
    for(name.comp in names.coef[-i])
    {
      p <- colnames(pca.loadings) == name.comp
      clp <- pca.loadings[l,p]
      bj <- pca.coef[names.coef == name.comp]
      s.j <- s.j + bj * clp / sdevs[names(m.s) == name.base]
    }
    b <- c(b, s.j)
    names(b)[length(b)] <- name.base
  }
  
  b
}

print(get.original.coef(x.red[,-1], pca.red, lm.red.pc1))

# Далі - прогнозуємо значення
# Переходимо до першого базису
x.test.full <- data.frame(predict(pca.full, x.dat.to.predict[,-1])[,1:P.full])
# Прогнозування
cl.predicted.full <- predict(lm.full, x.test.full)
# Обчислення залишків прогнозу
u.full <- x.dat.to.predict$cl - cl.predicted.full
# Переходимо до другого базису
x.test.red <- data.frame(predict(pca.red, x.dat.to.predict[,-1])[,1:P.red])
# Прогнозування
cl.predicted.red <- predict(lm.red, x.test.red)
cl.predicted.red.pc1 <- predict(lm.red.pc1, x.test.red)
cl.predicted.red.pc13 <- predict(lm.red.pc13, x.test.red)
# Обчислення залишків прогнозу
u.red <- x.dat.to.predict$cl - cl.predicted.red
u.red.pc1 <- x.dat.to.predict$cl - cl.predicted.red.pc1
u.red.pc13 <- x.dat.to.predict$cl - cl.predicted.red.pc13

# Однак перш ніж робити графік, проробимо кроки для прогнозування результатів
# з першої лабораторної роботи
times <- rownames(x.dat.to.predict)
rownames(x.dat.to.predict) <- NULL
lm.best.lab1 <- lm(cl ~ clx + cme - 1, data=x.red[-c(19,34,35),])
cl.predicted.lab1 <- predict(lm.best.lab1, x.dat.to.predict[,-1])
u.lab1 <- x.dat.to.predict$cl - cl.predicted.lab1

residuals.total <- c(u.full, u.red, u.red.pc1, u.red.pc13, u.lab1)
plot(times, u.full, ylim=c(min(residuals.total), max(residuals.total)+1), 
     type='l', col='darkgreen', main="Динаміка залишків на нових даних",
     xlab="Номер сесії")
lines(times, u.red, col='blue')
lines(times, u.red.pc1, col='magenta')
lines(times, u222, col="purple")
lines(times, u.red.pc13, col='orange')
lines(times, u.lab1, col="brown")
legend("topleft", 
       legend = list("PCA - повні дані", 
                     "PCA - 50 сесій - перші чотири компоненти", 
                     "PCA - 50 сесій - лише перша компонента",
                     "PCA - 50 сесій - лише перша та вилучено викид",
                     "PCA - 50 сесій - перша і третя компоненти",
                     "Lab1 - Найкраща модель"),
       col = c("darkgreen", "blue", "magenta", "purple", "orange", "brown"), lwd = 2)
abline(h=0, col='red', lty=2)
grid()

plot(lm.full)
qq.norm.intervals(MASS::studres(lm.full))

plot(lm.red.pc13)
qq.norm.intervals(MASS::studres(lm.red.pc13))