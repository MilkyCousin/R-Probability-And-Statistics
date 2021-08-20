library(nleqslv)

global.seed <- 1

rgaussmixt <- function(n, mu.1, mu.2, s.1, s.2, p)
{
  ind <- sample(c(1,2), n, replace = T, prob = c(p, 1 - p))
  u.sample <- rnorm(n, mean = c(mu.1, mu.2)[ind], sd = c(s.1, s.2)[ind])
}

est.mm <- function(x, m.1, m.2, s.2, p)
{
  est.r <- mean(x^2)/p - (1-p)/p * (s.2^2 + m.2^2) - m.1^2
  est.r
}

est.mm.0 <- function(x, m.1, m.2, s.2, p, x.0 = 1)
{
  mx.2 <- mean(x^2)
  mm.eq <- function(t)
  {
    sqrt(t) - sqrt(mx.2/p - (1-p)/p * (s.2^2 + m.2^2) - m.1^2)
  }
  r <- nleqslv(
    x.0, mm.eq, method = "Newton"
  )
  tt <- seq(-10, 10, 0.01)
  plot(tt, mm.eq(tt), type='l')
  abline(v = r$x)
  r$x
}

d.gaussmixt <- function(s.1, s.2, m.1, m.2, p, x)
{
  p*dnorm(x, m.1, s.1) + (1-p)*dnorm(x, m.2, s.2)
}

n <- 5000
u <- rgaussmixt(n, 1, 0, 0.05, 0.75, 0.6)
print(est.mm.0(u, 1, 0, 0.75, 0.6, x.0 = 1))
print(est.mm(u, 1, 0, 0.75, 0.6))
uu <- replicate(1000, est.mm.0(rgaussmixt(n, 1, 0, 0.05, 0.75, 0.6), 1, 0, 0.75, 0.6))