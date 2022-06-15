f <- function(t)
{
  t + sin(2 * t)
}

generate.x <- function(time, n, x0, theta, sigma)
{
  delta <- time / 2^n
  times <- (0:(2^n)) * delta
  increments.w <- rnorm(2^n, mean = 0, sd = sqrt(delta))
  increments.f <- f(times[-1]) - f(times[-(2^n - 1)])
  cumsum(c(x0, theta * increments.f + sigma * increments.w))
}

theta.est.last <- function(time, X.t)
{
  X.t[length(X.t)] / f(time)
}

f.prime <- function(t)
{
  1 + 2 * cos(2 * t)
}

theta.est.ls <- function(time, n, X.t)
{
  times <- (0:(2^n - 1)) * time / 2^n
  nominator <- sum(f.prime(times) * diff(X.t))
  denominator <- 3 * time + 2 * sin(2 * time) + sin(4 * time) / 2
  nominator / denominator
}

sigma.sq.est <- function(time, X.t)
{
  sum(diff(X.t)^2) / time
}

set.seed(0)

T.val <- 200
N <- 16

x0 <- 3
theta <- 4
sigma <- 3

B <- 1000
sigma.sq.val <- c()
for(b in 1:B)
{
  X.traj <- generate.x(T.val, N, x0, theta, sigma)
  sigma.sq.val <- c(sigma.sq.val, sigma.sq.est(T.val, X.traj))
}

