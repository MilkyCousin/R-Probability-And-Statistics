EulerSolverSDE <- function(T.val, n, a, b, x0)
{
  # Визначимо відставнь між вузлами
  delta <- T.val / n
  # Визначимо вузли
  t.nodes <- (0:n) * delta
  # Наближення у початковий момент часу
  X.approx <- x0
  # Рекурентне обчислення значень у наступних вузлах
  W.incr <- rnorm(n, mean = 0, sd = sqrt(delta))
  for(j in (1:length(t.nodes))[-1])
  {
    a.prev <- a(t.nodes[j-1], X.approx[j-1])
    b.prev <- b(t.nodes[j-1], X.approx[j-1])
    X.next <- X.approx[j-1] + a.prev * delta + b.prev * W.incr[j-1]
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
  1 + sin(2 * x)
}

b <- function(t, x)
{
  2 - cos(x)
}

theta <- -1

n <- 1000
x0 <- 1
B <- 1000

T.val <- 200

# Знаходимо 1000 різних траєкторій
X.times <- (0:n) * T.val / n
X.values.approx <- data.frame()

# Для повторного відтворення траєкторій явного розв'яку
set.seed(0)
for(k in 1:B)
{
  approx.result <- EulerSolverSDE(T.val, n, 
                                  function(t,x) {theta * a(t,x)}, b, x0)
  X.values.approx[1:(n+1),ncol(X.values.approx) + 1] <- approx.result$values
}

# Обчислення оцінки для зсуву
nominator <- function(X.traj, a, b)
{
  f <- function(x) { a(1, x) / (b(1, x))^2 }
  sum(f(X.traj[-length(X.traj)]) * diff(X.traj))
}

denominator <- function(X.traj, a, b, delta)
{
  f <- function(x) { (a(1, x) / b(1, x))^2 }
  sum(f(X.traj[-length(X.traj)]) * delta)
}

theta.mle <- function(X.traj, a, b, delta)
{
  nominator(X.traj, a, b) / denominator(X.traj, a, b, delta)
}

theta.estim <- apply(X.values.approx, 2, 
                     function(x) { theta.mle(x, a, b, T.val/n) })
names(theta.estim) <- NULL
print(c(mean(theta.estim), var(theta.estim)))