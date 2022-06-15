# Binom(p, 2); n = 65
# H0: p = 0.5
# H1: p = 0.6

set.seed(0)

inf.f <- function(m, p)
{
  m/(p*(1 - p))
}

ll.stat <- function(x, m, p.0, p.1)
{
  n <- length(x)
  sum(x)*log((p.1*(1 - p.0))/(p.0*(1 - p.1))) + m*n*log((1 - p.1)/(1 - p.0))
}

m <- 2
p.0 <- 0.5
p.1 <- 0.6 # 0.5 + 0.1 = 0.5 + v/sqrt(65)

n <- 65
v <- sqrt(n)*0.1

alpha <- 0.05

C.alpha <- -v^2/2*inf.f(m, p.0) + qnorm(1 - alpha) * v * sqrt(inf.f(m, p.0))
beta.pi <- pnorm((C.alpha - v^2/2*inf.f(m, p.0))/(v*sqrt(inf.f(m, p.0))))

alpha.0 <- alpha
beta.0 <- alpha
n.min <- (qnorm(1 - beta.0) + qnorm(1 - alpha.0))^2/(inf.f(m, p.0)*(p.1 - p.0)^2)

approx.ll.test <- function(x, m, p.0, p.1)
{
  ll.stat.x <- ll.stat(x, m, p.0, p.1)
  (ll.stat.x > C.alpha) * 1
}

B <- 10^4
ll.h.0 <- replicate(B, ll.stat(rbinom(n, m, p.0), m, p.0, p.1))
print(quantile(ll.h.0, 1 - alpha))
hist(ll.h.0)
