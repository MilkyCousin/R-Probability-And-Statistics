library("nleqslv")
set.seed(0)

w <- function(t)
{
  1 / (1 + t)
}

lambda <- 0.5

F.theor <- function(t)
{
  pexp(t, rate = lambda)
}

unbiased.sampling <- function(n)
{
  rexp(n, rate = lambda)
}

biased.sampling1 <- function(n, k, w)
{
  to.return <- c()
  while(length(to.return) != n)
  {
    generated <- rexp(1, rate = lambda)
    p <- k * w(generated)
    if(sample(c(T, F), 1, prob = c(p, 1 - p)))
    {
      to.return <- c(to.return, generated)
    }
  }
  to.return
}

elem <- function() {
  repeat {
    x <- rexp(1, lambda)
    if(runif(1) < w(x)) break
  }
  x
}

biased.sampling <- function(n, k, w)
{
  u <- replicate(n, elem())
  #print("finish")
  u
}

Femp <- function(x)
{
  f <- function(t) { sum(x < t) / length(x) }; f
}

FHT <- function(x, w)
{
  w.x <- w(x)
  f <- function(t) { sum(as.numeric(x < t) / w.x) / sum(1 / w.x) }; f
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
  f <- function(t) { sum(p * (z < t)) }; f
}

abs.dev <- function(F.theor, F.estim, nodes)
{
  max(abs(F.theor(nodes) - F.estim(nodes)))
}

n <- 300
x.unbiased <- unbiased.sampling(n)
x.biased <- biased.sampling(n, 1, w)

q1 <- qexp(0.01, lambda)
q2 <- qexp(0.99, lambda)
B <- 1000
I <- q2 * (0:B)/B + q1 * (B:0)/B

Fx.emp <- Femp(x.unbiased)
Fx.HT <- FHT(x.biased, w)
Fx.V <- FV(x.unbiased, x.biased, w, plot.c = T)

plot(I, sapply(I, Fx.emp), type = "s", col = "blue", 
     lty = 2, xlab = "t", ylab = "F(t)", main = "Емпірична функція розподілу")
lines(I, F.theor(I), col = "red", type = "l"); grid()
legend("bottomright", legend = c("Теоретична ф.р.", "Оцінка ф.р."),
       col = c("red", "blue"), lty = c(1, 2), lwd = 1)
print(paste("Empirical:", abs.dev(function(t) {sapply(t, Fx.emp)}, 
                                  function(t) {F.theor(t)}, I)))

plot(I, sapply(I, Fx.HT), type = "s", col = "blue", 
     lty = 2, xlab = "t", ylab = "F(t)", main = "Оцінка Горвіца-Томпсон")
lines(I, F.theor(I), col = "red", type = "l"); grid()
legend("bottomright", legend = c("Теоретична ф.р.", "Оцінка ф.р."),
       col = c("red", "blue"), lty = c(1, 2), lwd = 1)
print(paste("Horvitz-Thompson:", abs.dev(function(t) {sapply(t, Fx.HT)}, 
                                         function(t) {F.theor(t)}, I)))

plot(I, sapply(I, Fx.V), type = "s", col = "blue", 
     lty = 2, xlab = "t", ylab = "F(t)", main = "Оцінка Варді")
lines(I, F.theor(I), col = "red", type = "l"); grid()
legend("bottomright", legend = c("Теоретична ф.р.", "Оцінка ф.р."),
       col = c("red", "blue"), lty = c(1, 2), lwd = 1)
print(paste("Vardi:", abs.dev(function(t) {sapply(t, Fx.V)}, 
                              function(t) {F.theor(t)}, I)))

convF <- function(t) { 1/2 * (Fx.emp(t) + Fx.HT(t)) }
plot(I, sapply(I, convF), type = "s", col = "blue", 
     lty = 2, xlab = "t", ylab = "F(t)", main = "Опукла комбінація, lambda = 0.5")
lines(I, F.theor(I), col = "red", type = "l"); grid()
legend("bottomright", legend = c("Теоретична ф.р.", "Оцінка ф.р."),
       col = c("red", "blue"), lty = c(1, 2), lwd = 1)
print(paste("Convex:", abs.dev(function(t) {sapply(t, convF)}, 
                              function(t) {F.theor(t)}, I)))

