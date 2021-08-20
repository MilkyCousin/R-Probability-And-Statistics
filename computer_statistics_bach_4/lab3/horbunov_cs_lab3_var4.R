EstMM <- function(x, s, p)
{
  sqrt(mean(x^2)/p - s^2 * (1 - p)/p - 1)
}

seed.value <- 0

generate.mm <- function(theta, n = 1000, B = 1000, s = 3/4, p = 0.6)
{
  set.seed(seed.value)
  
  m.n <- c(1, 0)
  s.n <- c(theta, s)
  b.mm <- numeric(B)
  for(i in 1:B)
  {
    ind <- sample(c(1,2), n, prob = c(p, 1 - p), replace = T)
    b <- rnorm(n, mean = m.n[ind], sd = s.n[ind])
    b.mm[i] <- EstMM(b, s, p)
  }
  list(
    mm.est = b,
    asm.var = asymptot.var(s, p, theta),
    emp.var = n * var(b.mm),
    emp.bias = sqrt(n) * (mean(b) - theta)
  )
}

asymptot.var <- function(s, p, theta)
{
  a <- 3*p*theta^4+3*(1-p)*s^4+6*p^2*theta^2+6*p*(1-p)*s^2+p
  b <- 4 * p^2 * theta^2
  a/b
}

proba.get <- function(theta, s, p, alpha = 0.05, n = 2000, B = 1000)
{
  set.seed(seed.value)
  
  z <- qnorm(1 - alpha/2)
  
  frequencies <- replicate(
    B, {
      ind <- sample(c(1,2), n, prob = c(p, 1 - p), replace = T)
      b <- rnorm(n, mean = c(1, 0)[ind], sd = c(theta, s)[ind])
      mm.estimate <- EstMM(b, s, p)
      ifelse(abs(mm.estimate - theta) < sqrt(asymptot.var(s, p, theta)/n), 1, 0)
    }
  )
  1 - mean(frequencies)
}

for(theta.true in c(0.5, 1, 2, 4))
{
  results <- generate.mm(theta = theta.true)
  hist(
    results$mm.est, 
    probability = T, 
    main = paste("Histogram of estimates, sigma = ", theta.true),
    breaks = 20,
    xlab = "Value", ylim = c(0, .25 + dnorm(0))
  )
  curve(
    dnorm(
      x, mean = mean(results$mm.est), 
      sd = sd(results$mm.est)
      )
    , add = T, col = "red"
    )
}