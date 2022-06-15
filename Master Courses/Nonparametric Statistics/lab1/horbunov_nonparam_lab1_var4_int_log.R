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

FKM.var.estim <- function(t, x, d, FKM.at.t)
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
  (1 - FKM.at.t)^2 * sum(sums)
}

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

FKM.ci <- function(t, x, d, alpha)
{
  z <- qnorm(1 - alpha / 2)
  FKM.at.t <- FKM(t, x, d)
  sqrt.v <- sqrt(FKM.var.estim(t, x, d, FKM.at.t))
  FKM.at.t + z * sqrt.v * c(-1,1)
}

gFKM.ci <- function(t, x, d, alpha)
{
  z <- qnorm(1 - alpha / 2)
  FKM.at.t <- FKM(t, x, d)
  sqrt.s <- sqrt(gFKM.var.estim(t, x, d, FKM.at.t))
  1 - (1 - FKM.at.t) * exp(c(1,-1) * z * sqrt.s)
}

set.seed(0)

n <- 10

p <- 0.5
Qp <- qlnorm(p, meanlog = m.log, sdlog = s.log)

alpha <- 0.05

m <- 1000

fkm.ci.m <- c()
gfkm.ci.m <- c()

for(b in 1:m)
{
  u.boot <- gencens(n)
  Y.boot <- u.boot$dat
  d.boot <- u.boot$ind
  fkm.ci.b <- FKM.ci(Qp, Y.boot, d.boot, alpha)
  gfkm.ci.b <- gFKM.ci(Qp, Y.boot, d.boot, alpha)
  fkm.ci.m <- c(fkm.ci.m, (fkm.ci.b[1] > p) | (p > fkm.ci.b[2]))
  gfkm.ci.m <- c(gfkm.ci.m, (gfkm.ci.b[1] > p) | (p > gfkm.ci.b[2]))
}

alpha.fkm <- mean(fkm.ci.m)
alpha.gfkm <- mean(gfkm.ci.m)
print(c(alpha.fkm, alpha.gfkm))