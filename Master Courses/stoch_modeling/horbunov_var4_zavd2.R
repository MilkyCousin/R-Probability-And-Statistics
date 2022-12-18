# Моделювання траєкторії наближеного розв'язку СДР
# dX(t) = a(t,X(t))dt + b(t,X(t))dW(t), X(0) = x0, t in [0,T]
# Використовуючи метод Ейлера. Аргументи:
# x0 -- початкова точка розв'язку X(t)
# t -- T, верхній горизонт
# n -- кількість ненульових вузлів на [0,T]
# a, b -- коефіцієнти зсуву та дифузії відповідно
model.euler <- function(x0, t, n, a, b)
{
  # Ініціалізація кроку
  delta.n <- t/n
  x.approx <- numeric(length=n+1)
  # Старт наближення -- це старт розв'язку
  x.approx[1] <- x0
  # Буквально за означенням методу Ейлера
  for(j in 1:n)
  {
    dw <- rnorm(1, mean=0, sd=sqrt(delta.n))
    tk <- (j-1)*delta.n
    a.dt <- a(tk, x.approx[j]) * delta.n
    b.dw <- b(tk, x.approx[j]) * dw
    x.approx[j+1] <- x.approx[j] + a.dt + b.dw
  }
  x.approx
}

# Моделювання траєкторії розв'язку СДР
# dX(t) = 6X(t)dt + 1.6dW(t), X(0) = -3.8
model.solution <- function(t, n)
{
  # Ініціалізація кроку
  delta.n <- t/n
  # Обчислення рівномірної сітки
  time.set.n <- (0:n) * delta.n
  # Стартове значення розв'язку
  x0 <- -3.8
  # Обчислення експоненти
  exp.plus <- exp(6 * time.set.n)
  # Обчислення приростів Вінерівського процесу
  dw <- rnorm(n, mean=0, sd=sqrt(delta.n))
  # Наближене обчислення траєкторії інтеграла Іто
  exp.minus <- 1 / exp.plus[-(n+1)]
  int.wiener <- cumsum(c(0, exp.minus * dw))
  # Побудова траєкторії розв'язку на сітці
  x.approx <- (1.6 * int.wiener - 3.8) * exp.plus
  x.approx
}

#####

set.seed(777)

k <- 10
delta <- 2^(-k)
t <- 20
n <- t * 2^k
print(paste("n =", n))

# Кількість траєкторій наближеного розв'язку СДР для моделювання
n.copies <- 1000

# Безпосередньо моделювання траєкторій
x.copies <- matrix(nrow=n+1, ncol=n.copies)
for(k in 1:n.copies)
{
  x.traj <- model.solution(t, n)
  x.copies[,k] <- x.traj
}

#
par(mfrow=c(1,2))

time.set <- (0:n) * t/n
matplot(
  x=time.set, y=x.copies, type="l",
  main=paste("Real solution; n =", n, " delta.n =", round(t/n, 5))
)
grid()

# Безпосередньо моделювання траєкторій

y.copies <- matrix(nrow=n+1, ncol=n.copies)
for(k in 1:n.copies)
{
  y.traj <- model.euler(x0=-3.8, t, n, 
                        function(t,x) {6*x}, function(t,x) {1.6})
  y.copies[,k] <- y.traj
}

time.set <- (0:n) * t/n
matplot(
  x=time.set, y=y.copies, type="l",
  main=paste("Approximate solution; n =", n, " delta.n =", round(t/n, 5))
)
grid()

par(mfrow=c(1,1))

# Обчислення абсолютної похибки

abs.dev <- abs(x.copies - y.copies)

sq.dev.mean <- sapply(1:(n+1), function(j) {
  mean(abs.dev[j,]^2)
})

sq.dev.mean.sample <- sq.dev.mean[
  (time.set == 1) | (time.set == 10) | (time.set == 20)
]

n.specific <- 20 * 2^(3:5)
err.t1 <- c(14.23434, 13.69849, 12.96066)
err.t10 <- c(122.70524, 122.69308, 122.68198)
err.t20 <- c(242.70524, 242.69503, 242.69331)