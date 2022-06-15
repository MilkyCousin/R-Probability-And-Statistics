ClosedSolutionSDE <- function(times, w.incr)
{
  exponent <- cumsum(w.incr)
  x.closed <- 2 * c(1, exp(2 * exponent - times[-1])) 
  x.closed
}

EulerSolverSDE <- function(times, delta.noise, a, b, x0)
{
  delta <- diff(times)[1]
  n <- length(times) - 1
  x.approx <- x0
  for(j in (1:(n+1))[-1])
  {
    a.prev <- a(times[j-1], x.approx[j-1])
    b.prev <- b(times[j-1], x.approx[j-1])
    x.next <- x.approx[j-1] + a.prev * delta + b.prev * delta.noise[j-1]
    x.approx <- c(x.approx, x.next)
  }
  x.approx
}

# a(t,x) та b(t,x) з умови задачі
a <- function(t, x)
{
  x
}

b <- function(t, x)
{
  2*x
}

x0 <- 2
n <- 1000
B <- 1000

T.val <- 20

# Знаходимо 1000 різних траєкторій
delta <- T.val / n
X.times <- (0:n) * delta
X.values.approx <- data.frame()
X.values.closed <- data.frame()

# Для повторного відтворення траєкторій явного розв'яку
set.seed(10)
for(k in 1:B)
{
  w.incr <- rnorm(n, mean = 0, sd = sqrt(delta))
  approx.result <- EulerSolverSDE(X.times, w.incr, a, b, x0)
  closed.result <- ClosedSolutionSDE(X.times, w.incr)
  X.values.approx[1:(n+1),ncol(X.values.approx) + 1] <- approx.result
  X.values.closed[1:(n+1),ncol(X.values.closed) + 1] <- closed.result
}
rownames(X.values.approx) <- X.times
rownames(X.values.closed) <- X.times
matplot(X.times, X.values.approx, type="l", col=1:B, lty=1, ylim = c(0,max(X.values.approx)),
        main = paste("Траєкторії наближеного розв'язку при n =", n)); grid()
matplot(X.times, X.values.closed, type="l", col=1:B, lty=1, ylim = c(0,max(X.values.closed)),
        main = paste("Траєкторії явного розв'язку при n =", n)); grid()

interpolate <- function(t, times, X.values)
{
  n <- length(times)
  m <- times <= t
  j.nearest <- (1:n)[m][sum(m)]
  W.incr <- rnorm(1, mean = 0, sd = sqrt(t - time.idx))
  a.prev <- a(times[j.nearest], X.values[j.nearest])
  b.prev <- b(times[j.nearest], X.values[j.nearest])
  X.values[j.nearest] + a.prev * delta + b.prev * W.incr
}