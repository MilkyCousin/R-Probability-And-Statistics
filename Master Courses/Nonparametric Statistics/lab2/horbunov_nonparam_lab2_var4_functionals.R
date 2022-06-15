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

expected.value <- function(x, Fx, h = function(t) {t})
{
  sx <- sort(x)
  p <- diff(Fx(c(sx, Inf)))
  sum(p * h(sx))
}

sup.dist <- function()
{
  
}

set.seed(0)

# n <- 1000
# x.unbiased <- unbiased.sampling(n)
# x.biased <- biased.sampling(n, 1, w)
# 
# Fx.emp <- Femp(x.unbiased)
# Fx.HT <- FHT(x.biased, w)
# Fx.V <- FV(x.unbiased, x.biased, w)

# Математичне сподівання

# m.estim.emp <- expected.value(x.unbiased, Fx.emp)
# print(paste("Оцінка мат. сподівання за Femp:", m.estim.emp))
# m.estim.HT <- expected.value(x.biased, Fx.HT)
# print(paste("Оцінка мат. сподівання за FHT:", m.estim.HT))
# m.estim.conv <- 0.5 * (m.estim.emp + m.estim.HT)
# print(paste("Оцінка мат. сподівання за Fconv:", m.estim.conv))
# m.estim.V <- expected.value(c(x.unbiased, x.biased), Fx.V)
# print(paste("Оцінка мат. сподівання за FV:", m.estim.V))
  
imit.model.expectation <- function(n, k, w, h = function(t) {t}, B = 1000)
{
  exp.estim.tabl <- data.frame(empirical = double(),
                               horvitz.thompson = double(),
                               convex = double(),
                               vardi = double())
  for(b in 1:B)
  {
    x.repl.0 <- unbiased.sampling(n)
    x.repl.b <- biased.sampling(n, k, w)
    
    Fx.repl.emp <- Femp(x.repl.0)
    Fx.repl.HT <- FHT(x.repl.b, w)
    Fx.repl.V <- FV(x.repl.0, x.repl.b, w)
    
    m.repl.emp <- expected.value(x.repl.0, Fx.repl.emp)
    m.repl.HT <- expected.value(x.repl.b, Fx.repl.HT)
    m.repl.conv <- 0.5 * (m.repl.emp + m.repl.HT)
    m.repl.V <- expected.value(c(x.repl.0, x.repl.b), Fx.repl.V)
    
    exp.estim.tabl[nrow(exp.estim.tabl) + 1,] <- c(
      m.repl.emp, m.repl.HT, m.repl.conv, m.repl.V
    )
  }
  exp.estim.tabl
}

Ex <- 1 / lambda

N <- 1000
Ex.estim <- imit.model.expectation(N, 1, w)
Ex.estim.norm <- apply(Ex.estim, 2, function(v) { sqrt(N) * (v - Ex) })

bias.est <- apply(Ex.estim.norm, 2, mean)
var.est <- apply(Ex.estim.norm, 2, var)
for(j in 1:4)
{
  hist(Ex.estim.norm[,j], 
       prob = T, xlab = "value", 
       main = paste("Method:", colnames(Ex.estim.norm)[j], "for", "n:", N))
  curve(dnorm(x, bias.est[j], sqrt(var.est[j])), col = "red", add = T)
  grid()
}

print(N)
print(rbind(bias.est, var.est))