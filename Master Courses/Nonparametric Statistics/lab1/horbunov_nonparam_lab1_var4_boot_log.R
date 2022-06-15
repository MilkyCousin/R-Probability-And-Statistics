Femp <- function(t, x)
{
  sum(x < t) / length(x)
}

Fconf <- function(t, x, alpha = 0.05)
{
  F.emp.at.t <- Femp(t, x)
  z.alpha <- qnorm(1 - alpha / 2)
  var.at.t <- F.emp.at.t * (1 - F.emp.at.t) / length(x)
  ci.at.t <- F.emp.at.t + z.alpha * sqrt(var.at.t) * c(-1,1)
  ci.at.t
}

FKM <- function(t, x, d)
{
  n <- length(x)
  x.v <- sort(x)
  d.v <- d[order(x)]
  m <- x.v <= t
  if(!sum(m))
  {
    return(0)
  }
  idx <- (1:n)[m]
  mults <- sapply(idx, function(j) {
    1 - d.v[j] / (n - j + 1)
  })
  1 - prod(mults)
}

g <- function(t) { log(1 - t) }
g.prime <- function(t) { 1 / (t - 1) }

gFKM.var.estim <- function(t, x, d, FKM.at.t)
{
  n <- length(x)
  x.v <- sort(x)
  d.v <- d[order(x)]
  m <- x.v <= t
  if(!sum(m))
  {
    return(0)
  }
  idx <- (1:n)[m]
  sums <- sapply(idx, function(j) {
    d.v[j] / ((n - j + 1) * (n - j + 1 - d.v[j]))
  })
  (1 - FKM.at.t)^2 * sum(sums) * (g.prime(FKM.at.t))^2
}

m.log <- 0
s.log <- 1
df.cens <- 3

F.theor <- function(t)
{
  plnorm(t, meanlog = m.log, sdlog = s.log)
}

gencens <- function(N)
{
  t.vect <- rlnorm(N, meanlog = m.log, sdlog = s.log)
  c.vect <- rchisq(N, df = df.cens)
  d <- t.vect < c.vect
  z <- ifelse(d, t.vect, c.vect)
  z <- c(z, Inf)
  d <- c(d, 1)
  list(dat = z, ind = d)
}

set.seed(0)

n <- 1000

p <- 0.5
Qp <- qlnorm(p, meanlog = m.log, sdlog = s.log)

# ? 
m <- 1000
gFKM.boot.cor <- replicate(m, {
  u.boot <- gencens(n)
  Y.boot <- u.boot$dat
  d.boot <- u.boot$ind
  FKM.at.t <- FKM(Qp, Y.boot, d.boot)
  s.est <- sqrt(gFKM.var.estim(Qp, Y.boot, d.boot, FKM.at.t))
  (g(FKM.at.t) - g(p)) / s.est
})

hist.n <- hist(gFKM.boot.cor)
mh <- m * diff(hist.n$breaks)[1]
f <- function(t) { dnorm(t, mean = 0, sd = 1) * mh }
plot(hist.n, density = 15, col = "green", ylim = c(0, max(f(0), hist.n$counts)),
     xlab = "Значення", ylab = "Абсолютні частоти",
     main = paste("Гістограма для оцінки ln(1 - F(t)), n =", n))
curve(f(x), col = "red", add = T)
grid()