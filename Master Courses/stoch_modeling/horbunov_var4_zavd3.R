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

# Підрахунок оцінки методу найменших квадратів
numerator.ls <- function(x.sample, a, b)
{
  f <- function(x) { a(1, x) }
  sum(f(x.sample[-length(x.sample)]) * diff(x.sample))
}

denominator.ls <- function(x.sample, a, b, delta)
{
  f <- function(x) { a(1, x)^2 }
  sum(f(x.sample[-length(x.sample)]) * delta)
}

estim.ls <- function(x.sample, a, b, delta)
{
  N <- numerator.ls(x.sample, a, b)
  D <- denominator.ls(x.sample, a, b, delta)
  N / D
}

# Підрахунок оцінки методу максимальної правдоподібності
numerator.ml <- function(x.sample, a, b)
{
  f <- function(x) { a(1, x) / b(1, x)^2 }
  sum(f(x.sample[-length(x.sample)]) * diff(x.sample))
}

denominator.ml <- function(x.sample, a, b, delta)
{
  f <- function(x) { (a(1, x) / b(1, x))^2 }
  sum(f(x.sample[-length(x.sample)]) * delta)
}

estim.ml <- function(x.sample, a, b, delta)
{
  N <- numerator.ml(x.sample, a, b)
  D <- denominator.ml(x.sample, a, b, delta)
  N / D
}

# Для відтворення результатів
set.seed(788)

# Параметри СДР та наближення
x0 <- -1
a <- function(t, x) 
{
  -1 - 2*x
}
b <- function(t, x) 
{
  1 - sin(x) / 2
}

n <- 100
t <- 100

# Кількість траєкторій наближеного розв'язку СДР для моделювання
n.copies <- 1000

theta.real <- 2
# Безпосередньо моделювання траєкторій
x.copies <- matrix(nrow=n+1, ncol=n.copies)
for(k in 1:n.copies)
{
  x.traj <- model.euler(x0, t, n, function(t,x){ theta.real * a(t,x) }, b)
  x.copies[,k] <- x.traj
}

delta.n <- t/n
time.set <- (0:n) * delta.n
#matplot(
#  x=time.set, y=x.copies, type="l",
#  main=paste("n =", n, " delta.n =", round(t/n, 4))
#)
#grid()

theta.ls <- numeric(n.copies)
theta.ml <- numeric(n.copies)
for(k in 1:n.copies)
{
  x.cur.copy <- x.copies[,k]
  x.theta.ls <- estim.ls(x.cur.copy, a, b, delta.n)
  theta.ls[k] <- x.theta.ls
  x.theta.ml <- estim.ml(x.cur.copy, a, b, delta.n)
  theta.ml[k] <- x.theta.ml
}
print(paste("T =", t,"; n =", n, "; n.copies =", n.copies))
print(paste("Least squares:"))
print(c(mean(theta.ls) - theta.real, var(theta.ls)))
print(paste("Maximum likelihood:"))
print(c(mean(theta.ml) - theta.real, var(theta.ml)))