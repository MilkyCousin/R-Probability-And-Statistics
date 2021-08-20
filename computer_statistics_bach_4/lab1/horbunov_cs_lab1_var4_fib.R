library(rgl)

# Для повторного відтворення результатів генераторів R
set.seed(0)

# Аналогічні операції у перших двох файлах. Код дублюється, 
# можна було б виконати певні спрощення.

M <- 18

a.2 <- 75831
m.2 <- 2^31
i <- 2^15

rand <- function()
{
  i <<- (a.2 * i) %% m.2
  i
}

I <- replicate(M, rand())

l <- 17
k <- 5
m <- 2^31

fib.generator <- function()
{
  K <- (I[M - k] + I[M - l])%%m
  I <<- c(I[-1], K)
  K/m
}

lognormal.box.muller <- function(m = 0, s = 1)
{
  repeat{
    u <- fib.generator()
    v <- fib.generator()
    if((u != 0) && (v != 0))
    {
      g <- cos(2 * pi * u) * sqrt(-2 * log(v))
      return(exp(m + s * g))
    }
  }
}

f.div.cg <- function(t)
{
  (1/(t*sqrt(2*pi)))*(9/4)*(1+t^2)*exp(-8*(log(t)-1)^2/9)
}

lognormal.pdf.method <- function(m = 0, s = 1)
{
  repeat{
    c <- tan(pi * (fib.generator() - 0.5))
    u <- fib.generator()
    if(c > 0)
    {
      if(u < f.div.cg(c)) return(c)
    }
  }
}

plot.ecdf.theo.lognormal <- function(x, n, m = 0, s = 1, t = "")
{
  xs <- sort(x)
  
  plot(xs, ((1:n)/n), type='s', lwd=2,
       xlim=c(xs[1],xs[length(xs)]), ylim=c(0,1),
       xlab="t", ylab="F(t)", main=t)
  lines(xs, plnorm(xs, meanlog = m, sdlog = s), lty=3, lwd=3, col='red')
}

fast.lnorm.visual <- function(x.lnorm, n, m = 0, s = 1, t = "")
{
  plot.ecdf.theo.lognormal(x.lnorm, n, m, s, t)
  
  dn.lnorm <- dlnorm(x.lnorm, 1, 0.75)
  
  hist(x.lnorm, density=40, breaks=20, prob=TRUE,
       xlab="", ylab="", ylim=c(0, max(dn.lnorm)), main=t)
  curve(dlnorm(x.lnorm, 1, 0.75), col="red", lwd=2, 
        xname = "x.lnorm", add=TRUE, yaxt="n")
}

main <- function()
{
  N <- 500
  x <- replicate(500, fib.generator())
  xs <- sort(x)
  
  t <- "Custom + Fib Generator:"
  
  plot(xs, ((1:N)/N), type='s', lwd=2,
       xlim=c(0,1), ylim=c(0,1),
       xlab="t", ylab="F(t)", main=t)
  lines(xs, punif(xs), lty=3, lwd=3, col='red')
  
  plot(1:N, x, main=t, cex=0.45, xlab = 'index', ylab = 'value')
  
  x1 <- x[1:(N-2)]
  x2 <- x[2:(N-1)]
  x3 <- x[3:(N-0)]
  plot3d(x1,x2,x3)
  
  # ...
  
  log.mean <- 1
  log.sdev <- 0.75
  
  N <- 1000
  
  t.bm <- "Box-Muller transform + exponentiation method (Fib + Custom)"
  x.bm <- replicate(N, lognormal.box.muller(log.mean, log.sdev))
  fast.lnorm.visual(x.bm, N, log.mean, log.sdev, t.bm)
  
  t.pm <- "PDF Method (Fib + Custom)"
  x.pm <- replicate(N, lognormal.pdf.method(log.mean, log.sdev))
  fast.lnorm.visual(x.pm, N, log.mean, log.sdev, t.pm)
}