library(MASS)

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

# Ridge-регресія: підготовка

cv.vs.lambda <- function(ridge.model, j, text="")
{
  plot(ridge.model$lambda, ridge.model$GCV, type='l',
       xlab="lambda", ylab="CV", 
       main=paste("Графік залежності CV від lambda.", text))
  abline(v=ridge.model$lambda[j], col='red')
  grid()
}

coef.vs.lambda <- function(coefficients, lambda.set, j, text="", del.b0=T)
{
  m <- length(coefficients)
  matplot(lambda.set, coefficients[,-as.numeric(del.b0)], type="l",
          col=1:m, lty=1:m, lwd=2,
          xlab="lambda", ylab="coefficients",
          main=paste(
            "Динаміка значень коефіцієнтів при збільшенні lambda.", text
          ))
  abline(v=lambda.set[j], col='red', lwd=2)
  grid()
  legend("topright", legend=colnames(coefficients)[-as.numeric(del.b0)], 
         col=1:m, lty=1:m, lwd=1.5)
  
}

lambda.min <- 0.001
lambda.max <- 5 #50
lambdaset <- seq(lambda.min, lambda.max, 0.001)

# Ridge-регресія з штрафуванням зсуву

# Ridge-регресія за повними даними

x.dat.1 <- cbind(rep(1, nrow(x.dat)), x.dat)
colnames(x.dat.1)[1] <- "intercept"

# Підгонка
lm.ridge.full.1 <- lm.ridge(cl ~ . -1, data = x.dat.1, lambda = lambdaset)
l0 <- which.min(lm.ridge.full.1$GCV)

# Коефіцієнти моделі за оптимальним lambda, не нормовані
coef.full.l0 <- coef(lm.ridge.full.1)[l0,]

print(paste(
  "Коефіцієнти ridge-регресії за повними даними, lambda =", 
  lm.ridge.full.1$lambda[l0],
  "Штраф на зсув"
))
print(coef.full.l0)

# Графік залежності CV від lambda
#cv.vs.lambda(lm.ridge.full.1, l0, "Повні дані. Штраф наявний")

# Динаміка значень коефіцієнтів при збільшенні lambda
coef.vs.lambda(coef(lm.ridge.full.1), lambdaset, l0, 
               "Повні дані. Штраф наявний")

# Ridge-регресія за свіжими даними

x.red.1 <- cbind(rep(1, nrow(x.red)), x.red)
colnames(x.red.1)[1] <- "intercept"

lm.ridge.red.1 <- lm.ridge(cl ~ . -1, data = x.red.1, lambda = lambdaset)

# Номер lambda, який вважається оптимальним в сенсі зменшення значення CV
k0 <- which.min(lm.ridge.red.1$GCV)

# Коефіцієнти моделі за оптимальним lambda, не нормовані
coef.red.k0 <- coef(lm.ridge.red.1)[k0,]

print(paste(
  "Коефіцієнти ridge-регресії за свіжими даними, lambda =", 
  lm.ridge.red.1$lambda[k0],
  "Штраф на зсув"
))
print(coef.red.k0)

# Графік залежності CV від lambda
#cv.vs.lambda(lm.ridge.red.1, k0, "Свіжі дані. Штраф наявний")

# Динаміка значень коефіцієнтів при збільшенні lambda
coef.vs.lambda(coef(lm.ridge.red.1), lambdaset, k0, 
               "Свіжі дані. Штраф наявний")

################

# Ridge-регресія без штрафування зсуву

# Ridge-регресія за повними даними

# Підгонка
lm.ridge.full <- lm.ridge(cl ~ ., data = x.dat, lambda = lambdaset)

# Номер lambda, який вважається оптимальним в сенсі зменшення значення CV
j0 <- which.min(lm.ridge.full$GCV)

# Коефіцієнти моделі за оптимальним lambda, не нормовані
coef.full.j0 <- coef(lm.ridge.full)[j0,]
names(coef.full.j0)[1] <- "intercept"

print(paste(
  "Коефіцієнти ridge-регресії за повними даними, lambda =", 
  lm.ridge.full$lambda[j0],
  "Штраф на зсув відсутній"
))
print(coef.full.j0)

# Графік залежності CV від lambda
#cv.vs.lambda(lm.ridge.full, j0, "Повні дані. Штраф відустній")

# Динаміка значень коефіцієнтів при збільшенні lambda
coef.vs.lambda(coef(lm.ridge.full), lambdaset, j0, 
               "Повні дані. Штраф відустній")

# Ridge-регресія за свіжими даними

# Підгонка
lm.ridge.red <- lm.ridge(cl ~ ., data = x.red, lambda = lambdaset)

# Номер lambda, який вважається оптимальним в сенсі зменшення значення CV
i0 <- which.min(lm.ridge.red$GCV)

# Коефіцієнти моделі за оптимальним lambda, не нормовані
coef.red.i0 <- coef(lm.ridge.red)[i0,]
names(coef.red.i0)[1] <- "intercept"

print(paste(
  "Коефіцієнти ridge-регресії за свіжими даними, lambda =", 
  lm.ridge.red$lambda[i0],
  "Штраф на зсув відсутній"
))
print(coef.red.i0)

# Графік залежності CV від lambda
#cv.vs.lambda(lm.ridge.red, i0, "Свіжі дані. Штраф відустній")

# Динаміка значень коефіцієнтів при збільшенні lambda
coef.vs.lambda(coef(lm.ridge.red), lambdaset, i0, 
               "Свіжі дані. Штраф відустній")

##############

# Прогнозування

times <- rownames(x.dat.to.predict)
rownames(x.dat.to.predict) <- NULL
Y <- x.dat.to.predict$cl

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

# З поточної роботи
x.dat.to.predict.1 <- data.matrix(cbind(rep(1, N.test), x.dat.to.predict[,-1]))

Yhat_1_1 <- x.dat.to.predict.1%*%data.matrix(coef.full.l0)
Yhat_1_2 <- x.dat.to.predict.1%*%data.matrix(coef.full.j0)

Yhat_2_1 <- x.dat.to.predict.1%*%data.matrix(coef.red.k0)
Yhat_2_2 <- x.dat.to.predict.1%*%data.matrix(coef.red.i0)

U_1_1 <- drop(Y - Yhat_1_1)
U_1_2 <- drop(Y - Yhat_1_2)
U_2_1 <- drop(Y - Yhat_2_1)
U_2_2 <- drop(Y - Yhat_2_2)


U <- rbind(U_1_1, U_1_2, U_2_1, U_2_2, U.lab1, U.lab2)
nu <- ncol(U)
matplot(times, t(U), col=1:nu, lty=1,
        main="Динаміка залишків на нових даних",
        xlab="Номер сесії", ylab="Залишок",
        type='l')
grid()
abline(h=0, col='red', lty=2)
legend("topleft", legend=rownames(U), col=1:nu, lwd=1)