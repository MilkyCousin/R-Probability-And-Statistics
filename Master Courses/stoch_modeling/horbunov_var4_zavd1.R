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

seed.value <- 7
# 2000, 7, 128935
# Для відтворення результатів
set.seed(seed.value)

# Параметри СДР та наближення
x0 <- 4.4
a <- function(t, x) { 4^t * log(1 + x*x) }
b <- function(t, x) { (t + x) / (1 + t*t + x*x) }

n <- 100
t <- 1

# Кількість траєкторій наближеного розв'язку СДР для моделювання
n.copies <- 10

# Безпосередньо моделювання траєкторій
x.copies <- matrix(nrow=n+1, ncol=n.copies)
for(k in 1:n.copies)
{
  x.traj <- model.euler(x0, t, n, a, b)
  x.copies[,k] <- x.traj
}

time.set <- (0:n) * t/n
matplot(
  x=time.set, y=x.copies, type="l",
  main=paste("n =", n, " delta.n =", round(t/n, 4))
)
grid()
