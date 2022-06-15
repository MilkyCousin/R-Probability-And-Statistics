library("nleqslv")
set.seed(0)

w <- function(t)
{
  (cos(t / 2))^2
}

lambda <- 1

F.theor <- function(t)
{
  pexp(t, rate = lambda)
}

unbiased.sampling <- function(n)
{
  rexp(n, rate = lambda)
}

biased.sampling.old <- function(n, k, w)
{
  to.return <- c()
  while(length(to.return) != n)
  {
    generated <- rexp(1, rate = lambda)
    p <- k * w(generated)
    print(p)
    if(sample(c(T, F), 1, prob = c(p, 1 - p)))
    {
      to.return <- c(to.return, generated)
    }
  }
  to.return
}

biased.sampling <- function(n, k, w)
{
  to.return <- c()
  while(length(to.return) != n)
  {
    d <- n - length(to.return)
    generated <- rexp(d, rate = lambda)
    p <- k * w(generated)
    m <- sapply(p, function(q) { sample(c(T,F), 1, prob = c(q, 1 - q)) })
    to.return <- c(to.return, generated)
  }
  to.return
}

Femp <- function(x)
{
  f <- function(t) { sum(x < t) / length(x) }
  g <- function(v) { sapply(v, f) }; g
}

FHT <- function(x, w)
{
  w.x <- w(x)
  f <- function(t) { sum(as.numeric(x < t) / w.x) / sum(1 / w.x) }
  g <- function(v) { sapply(v, f) }; g
}

fv.target <- function(lambda, wz, n2)
{
  sum(wz / (lambda + n2 * wz)) - 1
}

FV <- function(unbiased, biased, w, x0 = 0, plot.c = F)
{
  n2 <- length(biased)
  z <- c(unbiased, biased)
  w.z <- w(z)
  to.solve <- function(l) { fv.target(l, w.z, n2) }
  solution <- nleqslv(x0, to.solve)
  lambda.opt <- solution$x
  if(plot.c)
  {
    t.val <- seq(0, (5/4 * lambda.opt), 0.01)
    plot(t.val, sapply(t.val, to.solve), type = 'l',
         xlab = "lambda", ylab = "sum")
    abline(h = 0, col = "black", lty = 2)
    abline(v = lambda.opt, col = "red", lty = 2)
    grid()
  }
  W <- 1 / (sum(1 / (lambda.opt + n2 * w.z)))
  p <- W / (lambda.opt + n2 * w.z)
  f <- function(t) { sum(p * (z < t)) }
  g <- function(v) { sapply(v, f) }; g
}

lambda.ivar.optim <- function(x.0, x.b, Fxemp, FxHT, w)
{
  sx.0 <- sort(x.0)
  sx.b <- sort(x.b)
  wsx <- w(sx.b)
  p <- diff(Fxemp(c(sx.0, Inf)))
  p.HT <- diff(FxHT(c(sx.b, Inf)))
  R.hat <- sum(wsx * p.HT)
  integral.emp <- 0
  I <- 1:length(x.0)
  for(k in I)
  {
    for(i in I[I < k])
    {
      for(j in I[I >= k])
      {
        integral.emp <- integral.emp + p[i]*p[j]*p[k]
      }
    }
  }
  pk.sub.wk <- p.HT / wsx
  integral.ht <- 0
  I <- 1:length(x.b)
  for(i in I)
  {
    res <- sum(pk.sub.wk[I < i] * (sum(p.HT[I >= k]))^2)
    res <- sum(pk.sub.wk[I >= i] * (sum(p.HT[I < k]))^2)
    integral.ht <- integral.ht + res * p.HT[i]
  }
  integral.ht <- R.hat * integral.ht
  lambda.opt <- integral.ht / (integral.ht + integral.emp)
  lambda.opt
}

n <- 300
x.unbiased <- unbiased.sampling(n)
x.biased <- biased.sampling(n, 1, w)

Fx.emp <- Femp(x.unbiased)
Fx.HT <- FHT(x.biased, w)
Fx.V <- FV(x.unbiased, x.biased, w)

lambda.opt.current <- lambda.ivar.optim(x.unbiased, x.biased, Fx.emp, Fx.HT, w)

# Прискористи biased.sampling
# Прискорити lambda.ivar.optim

B <- 1000
lambdas.opt <- c()
for(b in 1:1)
{
  print(b)
  x.0.boot <- unbiased.sampling(n)
  x.b.boot <- biased.sampling(n, 1, w)
  Fx.emp.boot <- Femp(x.0.boot)
  Fx.HT.boot <- FHT(x.b.boot, w)
  l.opt <- lambda.ivar.optim(x.0.boot, x.b.boot, Fx.emp.boot, Fx.HT.boot, w)
  lambdas.opt <- c(lambdas.opt, l.opt)
}