# Двокомпонентна гауссова суміш

d.norm.mixt <- function(t, m1, m2, s1, s2, p)
{
  p * dnorm(t, m1, s1) + (1 - p) * dnorm(t, m2, s2)
}

p.norm.mixt <- function(t, m1, m2, s1, s2, p)
{
  p * pnorm(t, m1, s1) + (1 - p) * pnorm(t, m2, s2)
}

q.norm.mixt <- function(alpha, m1, m2, s1, s2, p, q0 = 1)
{
  P <- function(q) { p.norm.mixt(q, m1, m2, s1, s2, p) - alpha }
  nleqslv::nleqslv(q0, P)$x
}

r.norm.mixt <- function(n, m1, m2, s1, s2, p)
{
  ind <- sample(c(1, 2), n, replace = T, prob = c(p, 1 - p))
  rnorm(n, mean = c(m1, m2)[ind], sd = c(s1, s2)[ind])
}

# Ядро Єпанєчнікова та відповідні числові характеристики від нього

epan.kernel <- function(t)
{
  3/4 * (1 - t^2) * (abs(t) < 1)
}

ep.d.sq <- 3/5
ep.D <- 1/5

# Ядерна оцінка щільності

dens.estim <- function(x, h, K)
{
  n <- length(x)
  g <- function(t)
  {
    sum(K((t - x) / h)) / (n * h)
  }
  gv <- function(t)
  {
    sapply(t, g)
  }
}

# Обчислення параметра згладжування за правилом Сільвермана

h.silv.rule <- function(d.sq, D, n, sigma.est)
{
  sigma.est * ((8 * sqrt(pi) * d.sq) / (3 * D^2 * n))^(1/5)
}

h.silv.simple <- function(d.sq, D, x)
{
  h.silv.rule(d.sq, D, length(x), sd(x))
}

h.silv.improved <- function(d.sq, D, x)
{
  h.silv.rule(d.sq, D, length(x), min(sd(x), IQR(x)/1.34))
}

# Теоретично оптимальний параметр згладжування

h.theor <- function(n)
{
  phi <- 1/(8*sqrt(pi))*(3/2-5/(2*exp(1)))
  (15/(phi * n))^(1/5)
}

# Непараметричне обчислення параметра згладжування - прямий підхід

phi.estim.epan <- function(x, h.pilot)
{
  n <- length(x)
  idx <- 1:n
  sx <- sort(x)
  s <- n * h.pilot + sum(sapply(idx[-1], function(i) {
    deltai <- sx[idx < i] - sx[i] + 2 * h.pilot
    sum(deltai * (deltai > 0))
  }))
  s * 9/2 * 1/(n^2 * h.pilot^6)
}

h.nonparam.epan <- function(x, h.pilot)
{
  phi.estim <- phi.estim.epan(x, h.pilot)
  (ep.d.sq / (length(x) * ep.D^2 * phi.estim))^(1/5)
}

# Непараметричне обчислення параметра згладжування через чисельні методи

dens.estim.epan.secderiv.sq <- function(x, h)
{
  n <- length(x)
  g <- function(t)
  {
    sum(abs(t - x) < h)^2 * 1 / (n^2 * h^6) * 9 / 4
  }
  gv <- function(t)
  {
    sapply(t, g)
  }
}

phi.estim.epan.calc <- function(x, h.pilot)
{
  f.secderiv.sq <- dens.estim.epan.secderiv.sq(x, h.pilot)
  integrate(
    f.secderiv.sq, min(x) - h.pilot, max(x) + h.pilot,
    subdivisions = 10^4
    )$value
}

h.nonparam.epan.calc <- function(x, h.pilot)
{
  phi.estim <- phi.estim.epan.calc(x, h.pilot)
  (ep.d.sq / (length(x) * ep.D^2 * phi.estim))^(1/5)
}

# Обчислення параметра згладжування за технікою кросс-валідації

f.a <- function(t, a)
{
  ((1 + t)^a - (-1)^a) / a
}

epan.kernel.conv <- function(t)
{ 
  ifelse(abs(t) >= 2, 0, {
    z <- -abs(t)
    (1 - z^2) * f.a(z, 1) + 2 * z * (f.a(z, 2) - f.a(z, 4)) + (z^2 - 2) * f.a(z, 3) + f.a(z, 5)
  })
}

CV.h <- function(h, x)
{
  n <- length(x)
  idx <- 1:n
  double.sum <- sum(sapply(idx, function(j) {
    delta <- (x[idx < j] - x[j])/h
    A <- sum(epan.kernel.conv(delta)) * 9/16
    B <- sum(epan.kernel(delta))
    A / n - 2 * B / (n - 1)
  }))
  (ep.d.sq + 2 * double.sum) / (n * h)
}

h.crossvalid <- function(x, h.min, h.max)
{
  optimize(function(h) { CV.h(h, x) }, c(h.min, h.max))$minimum
}

# Обчислення параметра згладжування за технікою кросс-валідації + FFT

g <- function(t)
{
  ifelse(abs(t) > 0, (sin(t) - t * cos(t)) / t^3, 1/3)
}

get.weights.ff <- function(a, b, m, x)
{
  n <- length(x)
  M <- 2^m
  delta <- (b-a)/M
  t <- a + (0:M) * delta
  weights <- rep(0, M)
  for(k in 1:(M-1))
  {
    yesno <- (t[k] <= x) & (x < t[k+1])
    weights[k] <- weights[k] + sum(ifelse(yesno, (t[k+1] - x)/(n*delta^2), 0))
    weights[k+1] <- weights[k+1] + sum(ifelse(yesno, (x - t[k])/(n*delta^2), 0))
  }
  weights
}

phi.at.zero.ff <- function(h, a, b, m, x)
{
  n <- length(x)
  M <- 2^m
  M2 <- M%/%2
  delta <- (b-a)/M
  t <- a + (0:M) * delta
  # M --> M + 1
  weights <- rep(0, M+1)
  for(k in 1:M)
  {
    yesno <- (t[k] <= x) & (x < t[k+1])
    weights[k] <- weights[k] + sum(ifelse(yesno, (t[k+1] - x)/(n*delta^2), 0))
    weights[k+1] <- weights[k+1] + sum(ifelse(yesno, (x - t[k])/(n*delta^2), 0))
  }
  s <- 2 * pi * (1:M2 + 1) / (b - a)
  Y <- (b - a) / (sqrt(2 * pi) * M) * fft(weights)  #FFT
  sum.val <- sum(g(h*s) * (3/sqrt(2*pi) * g(h*s) - 2) * abs(Y[1:M2 + 1])^2)
  1/(2*pi)^(3/2) * ((1 - 2*sqrt(2*pi))/2*sqrt(2*pi) + 6*(b-a)/M * sum.val)
}

CV.mod.ff <- function(h, a, b, m, x)
{
  phi.at.zero(h, a, b, m, x) + 2 * epan.kernel(0) / (length(x) * h)
}

h.crossvalid.ff <- function(x, h.min, h.max, a, b, m)
{
  optimize(function(s) { CV.mod(s, a, b, m, x) }, c(h.min, h.max))$minimum
}

# L2 - відстань

l2.sq <- function(f, g, a = -Inf, b = +Inf, n.subdiv = 10^4)
{
  i <- integrate(function(t) { (f(t) - g(t))^2 }, 
                 a, b, subdivisions = n.subdiv, stop.on.error = F)
  i$value
}

plot.l2.sq <- function(f, K, num.gen, N, seed.val = 0, n.sub = 10^3,
                       h.min = 0.01, h.max = 2, h.n = 200, B = 100)
{
  set.seed(seed.val)
  h.vals <- ((h.n-1):0)/(h.n-1) * h.min + (0:(h.n-1))/(h.n-1) * h.max
  l2.val <- rep(0, h.n)
  for(b in 1:B)
  {
    x.boot <- r.gen(N)
    for(j in 1:h.n)
    {
      l2.val[j] <- l2.val[j] + l2.sq(f, dens.estim(x.boot, h.vals[j], K))
    }
  }
  l2.val <- l2.val / B
  plot(h.vals, l2.val, type = "l", 
       xlab = "smoothing", ylab = "l2-squared", 
       main = paste("distance plot, N =", N, ", B =", B))
  grid()
}

# Моделювання

# Зернини: 0, 543787, 1912346 - те що випадково пробив клавішами, то й беру

set.seed(0)

N <- 500

m1 <- -1
m2 <- 1
s1 <- 1
s2 <- s1
p <- 0.5

f <- function(t) { d.norm.mixt(t, m1, m2, s1, s2, p) }
r.gen <- function(n) { r.norm.mixt(n, m1, m2, s1, s2, p) }

y <- r.norm.mixt(N, m1, m2, s1, s2, p)

# Графіки

alpha <- 0.01
D <- 1000

Qlower <- q.norm.mixt(alpha, m1, m2, s1, s2, p)
Qupper <- q.norm.mixt(1 - alpha, m1, m2, s1, s2, p)

I <- (D:0)/D * Qlower + (0:D)/D * Qupper

plot.dens <- function(dens.est, title)
{
  matplot(I,
          cbind(d.norm.mixt(I, m1, m2, s1, s2, p), dens.est(I)),
          col = c("red", "blue"),
          type = "l", xlab = "value", ylab = "density",
          main = paste(title, "N =", N))
  grid()
  legend("topleft", 
         legend = c("true density", "estimate"),
         col = c("red", "blue"), lwd = 1, lty = c(1,2))
}

# Правило Сільвермана - застосування

h.silv.s <- h.silv.simple(ep.d.sq, ep.D, y)
h.silv.i <- h.silv.improved(ep.d.sq, ep.D, y)

dens.silv.s <- dens.estim(y, h.silv.s, epan.kernel)
dens.silv.i <- dens.estim(y, h.silv.i, epan.kernel)

plot.dens(dens.silv.s, "silv rule - simple")
plot.dens(dens.silv.i, "silv rule - improved")

# Непараметричний підхід - застосування

h.nonp <- h.nonparam.epan.calc(y, h.silv.s)

dens.nonp.silv.s <- dens.estim(y, h.nonp, epan.kernel)

plot.dens(dens.nonp.silv.s, "nonparam - silv rule - simple")

# Теоретично оптимальний параметр згладжування - застосування

h.theor.val <- h.theor(N)
dens.theor <- dens.estim(y, h.theor.val, epan.kernel)
plot.dens(dens.theor, "theoretical")

# Кросс-валідація - застосування

h.mn <- 0.01
h.mx <- 2
#H <- seq(h.mn, h.mx, 0.05)
#plot(H, sapply(H, function(s) CV.h(s, y)), type = "l",
#     xlab = "h", ylab = "CV", main = paste("CV, N =", N))
h.cv <- h.crossvalid(y, h.mn, h.mx)
dens.cv <- dens.estim(y, h.cv, epan.kernel)
plot.dens(dens.cv, "cross-validation")

# Графіки MISE

#plot.l2.sq(f, epan.kernel, r.gen, N, seed.val = 0, 
#           n.sub = 10^4, h.min = h.mn, h.max = h.mx, h.n = 100, B = 100)
#abline(v = h.silv.i, col = "red", lty = 2)
#abline(v = h.nonp, col = "green", lty = 2)
#abline(v = h.cv, col = "purple", lty = 2)
#abline(v = h.theor.val, col = "blue", lty = 2)
#legend("topright", legend = c("Silverman", "Nonparametric", "CV", "Optimal"),
#       col = c("red", "green", "purple", "blue"), lwd = 1, lty = 2)

true.phi.val <- (3/2 - 5/(2*exp(1))) / (8 * sqrt(pi))

aMISE <- function(H)
{
  1/4 * true.phi.val * ep.D^2 * H^4 + ep.d.sq / H
}

H.min <- 0.01
H.max <- 10
H.seq <- seq(H.min, H.max, 0.01)
plot(H.seq, aMISE(H.seq), type = "l", xlab = "H", ylab = "aMISE(H)",
     main = paste("aMISE plot, N =", N))
grid()
abline(v = h.silv.i * N^(1/5), col = "red", lty = 2)
abline(v = h.nonp * N^(1/5), col = "green", lty = 2)
abline(v = h.cv * N^(1/5), col = "purple", lty = 2)
abline(v = h.theor.val * N^(1/5), col = "blue", lty = 2)
legend("topright", legend = c("Silverman", "Nonparametric", "CV", "Optimal"),
       col = c("red", "green", "purple", "blue"), lwd = 1, lty = 2)