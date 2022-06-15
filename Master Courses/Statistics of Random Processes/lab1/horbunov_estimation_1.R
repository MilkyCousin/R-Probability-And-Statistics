# Для повторного відтворення результатів
set.seed(0)

# Функція, яка повертає n-ий член послідовності bn.
# n - вектор з номерів, за якими потрібно обчислити ..
# .. відповідний елемент послідовності.
bn <- function(n) 
{ 
  ifelse(n > 0, 5 * n + 1 / n , 0) 
}

# Функція, яка повертає реалізацію випадкової послідовності у моделі вигляду:
# X(n) = x0 + theta * b(n) + sigma * R(n), n = 1, 2, ...
# X(0) = x0, b(0) = 0, R(n) = eps(1) + ... + eps(n), R(0) = 0, {eps(j)} - i.i.d.
# n - кількість перших n елементів послідовності (від 1 до n)
# x0 - початкове значення послідовності (в нульовий момент часу)
# theta - параметр зсуву послідовності
# sigma - параметр дифузії послідовності
# bn - перші n членів "трендової" послідовності
# eps - перші n членів "випадкової" послідовності
generate.sequence <- function(x0, theta, sigma, bn, eps)
{
  xn <- x0 + theta * bn + sigma * cumsum(eps)
  xn
}

# Обчислює оцінку зсуву за останнім спостереженням
theta.est.last <- function(n, xn, bn, x0) 
{
  (xn[n] - x0) / bn[n]
}

# Обчислює оцінку зсуву за МНК
theta.est.ls <- function(n, xn, bn, x0) 
{
  delta.b <- diff(c(0, bn))
  delta.x <- diff(c(x0, xn))
  sum(delta.x * delta.b) / sum(delta.b^2)
}

# Обчислює оцінку зсуву та початкового значення за МНК
theta.est.ls.2d <- function(xn, bn) 
{
  delta.b <- diff(bn)
  delta.x <- diff(xn)
  theta.ls <- sum(delta.x * delta.b) / sum(delta.b^2)
  x0.est <- xn[1] - theta.ls * bn[1]
  c(theta.ls, x0.est)
}

# Обчислює оцінку квадрата дифузії
sigma.estimate <- function(n, xn, x0)
{
  delta.x <- diff(c(x0, xn))
  sum(delta.x^2) / n
}

# Моделювання та обчислення

# Справжні значення параметрів
x0 <- -2
theta <- 0.4
sigma <- 2.3

# Кількість реалізацій процесу
B <- 1000

# Вектор з середніх по вибірці зі значень оцінки
theta.est.last.mean <- data.frame(norm = double(), unif = double())
# Вектор з дисперсій по вибірці зі значень оцінки
theta.est.last.var <- data.frame(norm = double(), unif = double())

# Вектор з середніх по вибірці зі значень нормованої оцінки
theta.est.last.mean.asympt <- data.frame(norm = double(), unif = double())
# Вектор з дисперсій по вибірці зі значень нормованої оцінки
theta.est.last.var.asympt <- data.frame(norm = double(), unif = double())

# Кількості перших значень, що спостерігаються
sizes <- c(50, 100, 200, 500, 1000)
mx.sizes <- max(sizes)
# Відповідно перші max(sizes) елементів послідовності b(n)
bn.full <- bn(1:mx.sizes)
# Для кожної кількості з sizes робимо ..
for(n in sizes)
{
  theta.est.repl <- data.frame(norm = double(), unif = double())
  # Моделювання по B реалізацій процесу за кожним розподілом похибок ..
  # .. та наближене обчислення відповідних числових характеристик
  for(b in 1:B)
  {
    x.repl.norm <- generate.sequence(x0, theta, sigma, bn.full[1:n], rnorm(n, 0, 1))
    x.repl.unif <- generate.sequence(x0, theta, sigma, bn.full[1:n], runif(n,-1, 1))
    theta.est.repl[nrow(theta.est.repl)+1,] <- c(
      theta.est.ls.2d(x.repl.norm, bn.full[1:n])[2] - x.repl.norm[1], # theta - 1, x0 - 2
      theta.est.ls.2d(x.repl.unif, bn.full[1:n])[2] - x.repl.unif[1]
    )
  }
  theta.est.last.mean[nrow(theta.est.last.mean)+1,] <- c(
    mean(theta.est.repl$norm),
    mean(theta.est.repl$unif)
  )
  theta.est.last.var[nrow(theta.est.last.var)+1,] <- c(
    var(theta.est.repl$norm),
    var(theta.est.repl$unif)
  )
}