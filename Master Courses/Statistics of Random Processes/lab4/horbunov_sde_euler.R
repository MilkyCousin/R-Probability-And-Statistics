EulerSolverSDE <- function(T.val, n, a, b, x0)
{
  # Визначимо відстань між вузлами
  delta <- T.val / n
  # Визначимо вузли
  t.nodes <- (0:n) * delta
  # Наближення у початковий момент часу
  X.approx <- x0
  # Рекурентне обчислення значень у наступних вузлах
  for(j in (1:length(t.nodes))[-1])
  {
    a.prev <- a(t.nodes[j-1], X.approx[j-1])
    b.prev <- b(t.nodes[j-1], X.approx[j-1])
    W.incr <- rnorm(1, mean = 0, sd = sqrt(delta))
    X.next <- X.approx[j-1] + a.prev * delta + b.prev * W.incr
    X.approx <- c(X.approx, X.next)
  }
  # Повертає вузли та значення у наближеного розв'язку СДР у них
  list(
    times = t.nodes,
    values = X.approx
  )
}

# a(t,x) та b(t,x) з умови задачі
a <- function(t, x)
{
  2 * x + 3 * sqrt(1 + t^2)
}

b <- function(t, x)
{
  2 * t - x / (1 + abs(x))
}

# Для повторного відтворення результатів
set.seed(0)

n <- 100000
x0 <- 3.2
B <- 10

T.val <- 1

# Знаходимо 10 різних траєкторій
X.times <- c()
X.values <- c()
for(k in 1:B)
{
  approx.result <- EulerSolverSDE(T.val, n, a, b, x0)
  X.times <- cbind(X.times, approx.result$times)
  X.values <- cbind(X.values, approx.result$values)
}
# Малюємо траєкторії
matplot(X.times, X.values, xlab = "t", ylab = "X(t)",
     type = "l", col = 1, ylim = c(min(X.values), max(X.values)),
     main = paste("Траєкторії наближеного розв'яку X(t), delta =", 1/n))
