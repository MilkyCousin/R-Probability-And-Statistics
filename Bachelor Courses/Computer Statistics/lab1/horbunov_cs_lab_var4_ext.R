I <- 2^15

a.1 <- 7^5
c.1 <- 0
m.1 <- 2^31-1

rand.uniform.pm <- function()
{
  # Генератор псевдовипадкових послідовностей на [0,1] з параметрами моделі
  # Парка та Мілера
  I <<- (a.1 * I + c.1) %% m.1
  I/m.1
}

lognormal.box.muller <- function(m = 0, s = 1)
{
  # Перетворення двох незалежних випадкових величин з рівномірним розподілом у
  # логнормальну величину за допомогою перетворення Бокса-Мюллера та взяття
  # експоненційної функції
  repeat{
    u <- rand.uniform.pm()
    v <- rand.uniform.pm()
    if((u != 0) && (v != 0))
    {
      g <- cos(2 * pi * u) * sqrt(-2 * log(v))
      return(exp(m + s * g))
    }
  }
}

f.div.cg <- function(t)
{
  # f(x)/(c*g(x)), c > 0 : f(x) <= c * g(x) для всіх x
  (1/(t*sqrt(2*pi)))*(9/4)*(1+t^2)*exp(-8*(log(t)-1)^2/9)
}

lognormal.pdf.method <- function(m = 0, s = 1)
{
  # Отримання логнормальної величини з використанням методу проріджування
  repeat{
    c <- tan(pi * (rand.uniform.pm() - 0.5)) # Cauchy variable
    u <- rand.uniform.pm() # U[0,1] variable
    if(c > 0)
    {
      if(u < f.div.cg(c)) return(c)
    }
  }
}

plot.ecdf.theo.lognormal <- function(x, n, m = 0, s = 1, t = "")
{
  # Викладки аналогічні тим, що наведені у попередньому файлі
  xs <- sort(x)
  
  plot(xs, ((1:n)/n), type='s', lwd=2,
       xlim=c(xs[1],xs[length(xs)]), ylim=c(0,1),
       xlab="t", ylab="F(t)", main=t)
  lines(xs, plnorm(xs, meanlog = m, sdlog = s), lty=3, lwd=3, col='red')
}

fast.lnorm.visual <- function(x.lnorm, n, m = 0, s = 1, t)
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
  N <- 1000
  
  log.mean <- 1
  log.sdev <- 0.75
  
  t.bm <- "Box-Muller transform + exponentiation method"
  x.bm <- replicate(N, lognormal.box.muller(log.mean, log.sdev))
  fast.lnorm.visual(x.bm, N, log.mean, log.sdev, t.bm)
  
  t.pm <- "PDF Method"
  x.pm <- replicate(N, lognormal.pdf.method(log.mean, log.sdev))
  fast.lnorm.visual(x.pm, N, log.mean, log.sdev, t.pm)
  
  t.sr <- "Standard R implementation"
  x.sr <- rlnorm(N, log.mean, log.sdev)
  fast.lnorm.visual(x.sr, N, log.mean, log.sdev, t.sr)
}