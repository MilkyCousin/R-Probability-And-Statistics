library(MASS)
library(leaps)

path <- "/home/coolpenguin/R/compsta/lab2"

source(paste(path, "qqplotinterval.R", sep="/"))

cur.data <- read.csv(paste(path, "merged.csv", sep="/"),
                     row.names = "X")

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

# Хід роботи

####

times <- rownames(x.dat.to.predict)
rownames(x.dat.to.predict) <- NULL
Y <- x.dat.to.predict$cl

####

r <- nrow(x.dat) - 1

# Відбір регресорів, повні дані
regsubs.full <- regsubsets(cl ~ ., data = x.dat,
                           nvmax = r, nbest = 5)

result.full <- summary(regsubs.full)
inclusions.full <- result.full$which
print(inclusions.full[1:5,])

Cp.full <- result.full$cp
p.full <- apply(inclusions.full, 1, sum)
plot(p.full, Cp.full)
abline(0, 1, col="red")
grid()

plot(p.full, Cp.full, ylim=c(0, 30))
abline(0, 1, col="red")
grid()

print(inclusions.full[36,])

lm.full.36 <- lm(cl ~ . - cmcsa, data=x.dat)

U.full.36 <- Y - predict(lm.full.36, x.dat.to.predict[,-1])

# Відбір регресорів, свіжі дані
regsubs.red <- regsubsets(cl ~ ., data = x.red,
                           nvmax = r, nbest = 5)

result.red <- summary(regsubs.red)
inclusions.red <- result.red$which
print(inclusions.red[1:5,])

Cp.red <- result.red$cp
p.red <- apply(inclusions.red, 1, sum)
plot(p.red, Cp.red)
abline(0, 1, col="red")
grid()

plot(p.red, Cp.red, xlim=c(3, 10), ylim=c(3,15))
abline(0, 1, col="red")
grid()

print(inclusions.red[26,])
print(inclusions.red[31,])

lm.red.26 <- lm(cl ~ . - clf - cmg - cmi, data=x.red)
lm.red.31 <- lm(cl ~ . - clf - cmg, data=x.red)

U.red.26 <- Y - predict(lm.red.26, x.dat.to.predict[,-1])
U.red.31 <- Y - predict(lm.red.31, x.dat.to.predict[,-1])

# З першої роботи
lab1.best <- lm(cl ~ clx + cme - 1, data=x.red[-c(19,34,35),])
U.lab1 <- Y - predict(lab1.best, x.dat.to.predict[,-1])

# З другої роботи
pca.red <- princomp(x.red[,-1], cor = T)
P.red <- 4
x.red.pc <- data.frame(pca.red$scores[,1:P.red])
x.dat.to.predict.pc <- predict(pca.red, x.dat.to.predict[,-1])

lab2.best <- lm(x.red$cl ~ Comp.1, data=x.red.pc)
U.lab2 <- Y - predict(lab2.best, data.frame(x.dat.to.predict.pc))

# З третьої роботи
lambda.i0 <- 0.236
lab3.best <- lm.ridge(cl ~ ., data = x.red, lambda = lambda.i0)
coef.ridge <- coef(lab3.best)

x.dat.to.predict.1 <- data.matrix(cbind(rep(1, N.test), x.dat.to.predict[,-1]))
U.lab3 <- drop(Y - x.dat.to.predict.1%*%data.matrix(coef.ridge))

#######

U <- rbind(U.full.36, U.red.26, U.red.31, U.lab1, U.lab2, U.lab3)
nu <- ncol(U)
colnums <- 1:nu+8
matplot(times, t(U), col=colnums, lty=1,
        main="Динаміка залишків на нових даних",
        xlab="Номер сесії", ylab="Залишок",
        type='l')
grid()
abline(h=0, col='red', lty=2)
legend("topleft", legend=c(
  "Повні дані, вилучено змінну cmcsa",
  "Свіжі дані, вилучено змінні clf, cmg, cmi",
  "Свіжі дані, вилучено змінні clf, cmg",
  "Найкраща модель з першої роботи",
  "Найкраща модель з другої роботи",
  "Найкраща модель з третьої роботи"
), col=colnums, lwd=1)