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

gencens <- function(N)
{
  t.vect <- rlnorm(N, meanlog = 0, sdlog = 1)
  c.vect <- rchisq(N, df = 3)
  d <- t.vect < c.vect
  z <- ifelse(d, t.vect, c.vect)
  z <- c(z, Inf)
  d <- c(d, 1)
  list(dat = z, ind = d)
}

set.seed(0)

n <- 1000

m.log <- 0
s.log <- 1

u <- gencens(n)

Y <- u$dat
d <- u$ind

FKM.Y <- function(t) { FKM(t, Y, d) }

q1 <- qlnorm(0.01, meanlog = m.log, sdlog = s.log)
q2 <- qlnorm(0.99, meanlog = m.log, sdlog = s.log)
D <- 1000
I <- q2 * (0:D)/D + q1 * (D:0)/D

FKM.vals <- sapply(I, FKM.Y)
FLN.vals <- plnorm(I, meanlog = m.log, sdlog = s.log)

plot(I, FKM.vals, type="s", ylim = c(0,1),
     xlab = "t", ylab = "F(t)", col = "blue",
     main = paste("Графіки функцій розподілу, n =", n))
lines(I, FLN.vals, col="red", lty=2)
grid()
legend("bottomright", legend = c("Оцінка Каплана-Мейєра", "Теоретична ф.р."),
       col = c("blue", "red"), lty = c(1, 2), lwd = 1)
