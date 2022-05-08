library("nleqslv")

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

biased.sampling <- function(n, k, w)
{
  to.return <- c()
  while(length(to.return) != n)
  {
    d <- n - length(to.return)
    generated <- rexp(d, rate = lambda)
    p <- k * w(generated)
    m <- sapply(p, function(q) { sample(c(T,F), 1, prob = c(q, 1 - q)) })
    to.return <- c(to.return, generated[m])
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

sup.dist <- function(x, Fest, Ftheor)
{
  sorted.x <- sort(x)
  fj <- Ftheor(sorted.x)
  hat.fj <- Fest(sorted.x)
  fplus <- max(abs(fj - hat.fj))
  fplus
}

set.seed(0)

n <- 300
x.unbiased <- unbiased.sampling(n)
x.biased <- biased.sampling(n, 1, w)

Fx.emp <- Femp(x.unbiased)
Fx.HT <- FHT(x.biased, w)
Fx.V <- FV(x.unbiased, x.biased, w)
Fx.conv <- function(t) { 0.5 * (Fx.emp(t) + Fx.HT(t)) }

#print(sup.dist(x.unbiased, Fx.emp, F.theor))
#print(sup.dist(x.biased, Fx.HT, F.theor))
#print(sup.dist(c(x.unbiased, x.biased), Fx.conv, F.theor))
#print(sup.dist(c(x.unbiased, x.biased), Fx.V, F.theor))

imit.model.distance <- function(n, k, w, B = 1000)
{
  dist.estim.tabl <- data.frame(empirical = double(),
                                horvitz.thompson = double(),
                                convex = double(),
                                vardi = double())
  for(b in 1:B)
  {
    x.repl.0 <- unbiased.sampling(n)
    x.repl.b <- biased.sampling(n, k, w)
    
    Fx.repl.emp <- Femp(x.repl.0)
    Fx.repl.HT <- FHT(x.repl.b, w)
    Fx.repl.conv <- function(t) { 0.5 * (Fx.repl.emp(t) + Fx.repl.HT(t)) }
    Fx.repl.V <- FV(x.repl.0, x.repl.b, w)
    
    d.emp <- sup.dist(x.repl.0, Fx.repl.emp, F.theor)
    d.HT <- sup.dist(x.repl.b, Fx.repl.HT, F.theor)
    d.conv <- sup.dist(c(x.repl.0, x.repl.b), Fx.repl.conv, F.theor)
    d.V <- sup.dist(c(x.repl.0, x.repl.b), Fx.repl.V, F.theor)
    
    dist.estim.tabl[nrow(dist.estim.tabl) + 1,] <- c(d.emp, d.HT, d.conv, d.V)
  }
  dist.estim.tabl
}

N <- 250
dist.estim.tabl.N <- imit.model.distance(N, 1, w, B = 1000)
print(paste("N =", N))
print(apply(dist.estim.tabl.N, 2, mean))