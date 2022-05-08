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

set.seed(0)

n <- 1000

lambda <- 1

Y <- rexp(n, lambda)
FY <- function(t) { Femp(t,Y) }

q1 <- qexp(0.01, lambda)
q2 <- qexp(0.99, lambda)
D <- 1000
I <- q2 * (0:D)/D + q1 * (D:0)/D

abs.dev <- abs((1:n)/n - pexp(sort(Y), lambda))
print(paste("Відхилення:", max(abs.dev)))

plot(I, sapply(I, FY), type="s", 
     xlab = "t", ylab = "F(t)", col = "blue",
     main = paste("Графіки функцій розподілу, n =", n))
curve(pexp(x, lambda), col="red", lty=2, add=T)
grid()
legend("bottomright", legend = c("Емпірична ф.р.", "Теоретична ф.р."),
       col = c("blue", "red"), lty = c(1, 2), lwd = 1)

p <- 1/3
q <- 1 - p
Qp <- qexp(p, lambda)
Qq <- qexp(q, lambda)

alpha <- 0.05
m <- 1000
p.cnt <- c()
q.cnt <- c()

for(b in 1:m)
{
  Yb <- rexp(n, lambda)
  CIp <- Fconf(Qp, Yb, alpha)
  CIq <- Fconf(Qq, Yb, alpha)
  p.cnt <- c(p.cnt, (CIp[1] < p) & (p < CIp[2]))
  q.cnt <- c(q.cnt, (CIq[1] < q) & (q < CIq[2]))
}

print(paste("Для p:", 1 - mean(p.cnt)))
print(paste("Для q:", 1 - mean(q.cnt)))