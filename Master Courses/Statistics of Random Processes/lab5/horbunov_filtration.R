a <- 1.3
b <- 0.8
A <- -1
B <- 1.5

k <- A^2 / B^2
d <- sqrt(a^2 + k * b^2)

h1 <- (a + d)/k
h2 <- (a - d)/k

c <- 1 / (h1 - h2)

integr.e <- function(t) { h1 * t + B^2 * log((exp(-t/(B^2 * c)) - h2 / h1) / (1 - h2 / h1)) }

gu <- function(t) { exp(a * t - k * integr.e(t)) }

g <- function(t, s) { gu(t) / gu(s) * (s <= t) }

h <- function(t) { (1 - exp(-t/(B^2 * c))) / (1 / h1 - 1 / h2 * exp(-t/(B^2 * c))) }

# G <- function(t, s) { A / B^2 * h(s) * g(t, s) }

# Моделювання X1

gen.time.nodes <- function(T.val, n)
{
  delta <- T.val / n
  (0:n) * delta
}

# Далі припускається, що t{j+1} - t{j} = delta для довільних j = 0, ..., n-1

# Моделювання X1, X2

model.x1 <- function(ts)
{
  delta <- diff(ts)[1]
  n <- length(ts) - 1
  w1.incr <- rnorm(n, mean = 0, sd = sqrt(delta))
  c(0, b * exp(a * ts[-1]) * cumsum(exp(-a * ts[-(n+1)]) * w1.incr))
}

model.x2 <- function(ts, x1.traj)
{
  delta <- diff(ts)[1]
  n <- length(ts) - 1
  w2.incr <- rnorm(n, mean = 0, sd = sqrt(delta))
  c(0, A * cumsum(x1.traj[-(n+1)]) * delta + B * cumsum(w2.incr))
}

# Моделювання оптимального фільтру для X1

model.filter <- function(ts, x2.traj)
{
  n <- length(ts) - 1
  delta.x2 <- diff(x2.traj)
  c(0, A / B^2 * gu(ts[-1]) * cumsum(h(ts[-(n+1)]) * delta.x2 / gu(ts[-(n+1)])))
}

aa <- function(t, x) { (a - k * h(t)) * x }
bb <- function(t, x) { A / B^2 * h(t) }

x1.hat.approx <- function(ts, x2)
{
  delta <- diff(ts)[1]
  n <- length(ts) - 1
  X.approx <- 0
  x2.incr <- diff(x2)
  for(j in (1:(n+1))[-1])
  {
    a.prev <- aa(ts[j-1], X.approx[j-1])
    b.prev <- bb(ts[j-1], X.approx[j-1])
    X.next <- X.approx[j-1] + a.prev * delta + b.prev * x2.incr[j-1]
    X.approx <- c(X.approx, X.next)
  }
  # Повертає значення у наближеного розв'язку СДР у них
  X.approx
}

# Практичне застосування

# set.seed(0)

T.val <- 1
m <- 8
n <- 2^m
times <- gen.time.nodes(T.val, n)

X1 <- model.x1(times)
X2 <- model.x2(times, X1)

X1.hat <- model.filter(times, X2)
X1.hat.approx <- x1.hat.approx(times, X2)

y.min <- min(X1, X2, X1.hat, X1.hat.approx)
y.max <- max(X1, X2, X1.hat, X1.hat.approx)
plot(times, X2, type="l", ylim = c(y.min, y.max), col = "darkgreen",
     xlab = "t", ylab = "X(t)", lwd = 2)
lines(times, X1, type="l", col = "red")
lines(times, X1.hat, type="l", col = "darkblue", lty = 2)
lines(times, X1.hat.approx, type="l", col = "purple", lty = 2)
grid()
legend("bottomleft", col = c("darkgreen", "red", "darkblue", "purple"),
       legend = c("X2", "X1", "X1.hat", "X1.hat.approx"), 
       lwd = c(2,1,1,1), lty = c(1,1,2,2))
stop("Будівництво")
print("Моделювання sigma.sq.t")

plot(times, h(times), type = "l")
B <- 10^4
sq.dist.repl <- replicate(B, {
  X1.new <- model.x1(times)
  X2.new <- model.x2(times, X1.new)
  X1.hat.new <- model.filter(times, X2.new)
  X1.hat.approx.new <- x1.hat.approx(times, X2.new)
  c((X1.new - X1.hat.new)^2, (X1.new - X1.hat.approx.new)^2)
})
estim.sigma.sq.t.int <- apply(sq.dist.repl[1:(n+1),], 1, mean)
estim.sigma.sq.t.appr <- apply(sq.dist.repl[(n+2):(2*(n+1)),], 1, mean)
lines(times, estim.sigma.sq.t.int, col = "blue", lty = 2)
lines(times, estim.sigma.sq.t.appr, col = "purple", lty = 2)
grid()
legend("topleft", col = c("black", "blue", "purple"),
       legend = c("Справжнє значення", 
                  "Оцінка за прямим обчисленням", 
                  "Оцінка за Ейлеровим наближенням"), 
       lwd = 1, lty = c(1,2,2))
