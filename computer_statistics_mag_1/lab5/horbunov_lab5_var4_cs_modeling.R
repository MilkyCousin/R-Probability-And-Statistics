# Фіксуємо зернину
set.seed(1)
#set.seed(10000)
# Функція "кореня з дисперсії"
g <- function(t)
{
  sqrt(0.5)*abs(t-1)
}
# Параметри зсуву і нахилу відповідно
a <- -2
b <- 0.5
# Мат. сподівання, дисперсія X ~ N(m, s)
m <- -1
s <- sqrt(2)
# Обсяг вибірки
N <- 1000
# Генеруємо дані
# Регресор - набір з реалізацій X
X <- data.matrix(rnorm(N, mean=m, sd=s))
# Похибки
E <- data.matrix(rnorm(N, mean=0, sd=g(X)))
# Моделюємо відгук
Y <- a + b*X + E

# Що, власне, ми створили?
plot(X, Y, cex=0.75, main="Діаграма розсіювання Y ~ X", xlab="X", ylab="Y")
grid()

# Класичний МНК
mX = mean(X)
mY = mean(Y)

b_ols = cov(X, Y)/var(X)
a_ols = mY - b_ols * mX
b <- data.matrix(c(a_ols, b_ols))

abline(a_ols, b_ols, col="red")
print(c(a_ols, b_ols))

# Прогноз-залишки
Yhat_ols <- cbind(1, X)%*%b
U_ols <- Y - Yhat_ols
plot(Yhat_ols, U_ols, main="Прогноз-Залишки", xlab="Yhat", ylab="U")
abline(h=0, col="red")
grid()

# Діаграма розсіювання квадрату залишків відносно регресора
plot(X, U_ols^2, main="Діаграма розсіювання", xlab="X", ylab="U^2")
grid()

# Підгонка параметрів у моделі U^2 = alpha + beta*X + gamma*X^2
X.poly <- cbind(1,X,X^2)
A.poly <- t(X.poly)%*%X.poly
greeks <- solve(A.poly)%*%t(X.poly)%*%(U_ols^2)
print(greeks)

f <- function(t) {c(1, t, t^2)%*%greeks}
fv <- function(v) {sapply(v, f)}
curve(fv(x), col="red", add=T)

Usq_fitted <- X.poly%*%greeks
print(sum(Usq_fitted < 0))

# Лінії рівня
Jw <- function(b, W)
{
  b <- data.matrix(b)
  t(Y - cbind(1, X)%*%b)%*%W%*%(Y - cbind(1, X)%*%b)
}

W.res <- diag(as.numeric(1 / Usq_fitted))

Jw.res <- function(bx, by) 
{
  apply(cbind(bx, by), 1, function(t) {Jw(t, W.res)})
}

#Ix <- seq(-4, 4, 0.1)
#Iy <- seq(-4, 4, 0.1)
#Jz <- outer(Ix, Iy, Jw.res)
#contour(Ix, Iy, Jz, nlevels=45, xlab="x", ylab="y", 
#        main="Лінії рівня навантаженого функціоналу МНК")

Aw <- t(cbind(1,X))%*%W.res%*%cbind(1,X)
bw <- solve(Aw)%*%t(cbind(1,X))%*%W.res%*%Y
#points(bw[1], bw[2], col="green")

#legend("bottomright", legend = "Критична точка функціоналу", 
#       col = "green", lwd = 1)

# Корекція на ваги

#W.res.0 <- W.res * (W.res > 0)

#Jw.res.0 <- function(bx, by) 
#{
#  apply(cbind(bx, by), 1, function(t) {Jw(t, W.res.0)})
#}

#Jz0 <- outer(Ix, Iy, Jw.res.0)
#contour(Ix, Iy, Jz0, nlevels=45, xlab="x", ylab="y", 
#        main="Лінії рівня виправленого навантаженого функціоналу МНК")

#Aw0 <- t(cbind(1,X))%*%W.res.0%*%cbind(1,X)
#bw0 <- solve(Aw0)%*%t(cbind(1,X))%*%W.res.0%*%Y
#points(bw0[1], bw0[2], col="green")

#legend("bottomright", legend = "Критична точка функціоналу", 
#       col = "green", lwd = 1)

Xw <- W.res%*%(X - mX)
Yw <- W.res%*%(Y - mY)
plot(Xw, Yw, main="Діаграма розсіювання перетворених X та Y",
     xlab="Xw", ylab="Yw")
grid()
abline(a=0, bw[2], col="red")

plot(Xw*bw[2], Yw - Xw*bw[2], main="Прогноз-Залишки",
     xlab="Yw_hat", ylab="Uw")
grid()
abline(h=0, col="red")

# Нормовані вагові коефіцієнти

W.res.norm <- W.res*(W.res <= 1) + 1*(W.res > 1)
Awn <- t(cbind(1,X))%*%W.res.norm%*%cbind(1,X)
bwn <- solve(Awn)%*%t(cbind(1,X))%*%W.res.norm%*%Y

Xwn <- W.res.norm%*%(X - mX)
Ywn <- W.res.norm%*%(Y - mY)
plot(Xwn, Ywn, main="Діаграма розсіювання перетворених X та Y",
     xlab="Xw", ylab="Yw")
grid()
abline(a=0, bwn[2], col="red")

plot(Xwn*bwn[2], Ywn - Xwn*bwn[2], main="Прогноз-Залишки",
     xlab="Yw_hat", ylab="Uw")
grid()
abline(h=0, col="red")