euler.approx <- function(times, w.incr, delta, a, b, x0)
{
  x.approx <- c(x0)
  for(j in (1:length(times))[-1])
  {
    A <- a(times[j], x.approx[j-1]) * delta
    B <- b(times[j], x.approx[j-1]) * w.incr[j-1]
    x.next <- x.approx[j-1] + A + B
    x.approx <- c(x.approx, x.next)
  }
  x.approx
}

a <- function(t, x)
{
  2 * x + 3 * sqrt(1 + t^2)
}

b <- function(t, x)
{
  2 * t - x / (1 + abs(x))
}

###

set.seed(0)

n <- 1000
x0 <- 3.2

B <- 10
x.approxed <- c()

delta <- 1/n
T.val <- 1
times <- (0:n) * delta * T.val

for(u in 1:B)
{
  w.incr <- rnorm(n, mean = 0, sd = sqrt(delta))
  x.approxed <- cbind(x.approxed, euler.approx(times, w.incr, delta, a, b, x0))
}

matplot(times, (x.approxed), xlab = "t", ylab = "X(t)",
        type = "l", col = 1, ylim = c(min(x.approxed), max(x.approxed)),
        main = paste("Траєкторії наближеного розв'яку X(t), delta =", 1/n))