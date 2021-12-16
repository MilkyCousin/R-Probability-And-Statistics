qq.norm.intervals <- function(x, B=1000, alpha=0.05, apply.intervals=T)
{
  n <- length(x)
  rep.samples <- matrix(rnorm(n*B), nrow=n, ncol=B)
  rep.samples <- apply(rep.samples, 2, sort)
  upper.bound <- apply(rep.samples, 1, quantile, probs=1-alpha/2)
  lower.bound <- apply(rep.samples, 1, quantile, probs=alpha/2)
  
  empirical.quantiles <- sort(x)
  theoretical.quantiles <- qnorm((1:n-0.5)/n)
  
  plot(rep(theoretical.quantiles, 3), 
       c(lower.bound, empirical.quantiles, upper.bound), type="n",
       xlab = "Theoretical quantiles", ylab = "Empirical quantiles",
       main = paste("QQ-diagram, std normal distrib, alpha = ", alpha))
  if(apply.intervals)
  {
    segments(x0=theoretical.quantiles, y0=lower.bound,
             x1=theoretical.quantiles, y1=upper.bound, col="purple")
  }
  abline(0, 1, lty=2, lwd=1 ,col='red')
  points(theoretical.quantiles, empirical.quantiles, cex=1, lwd=1)
  grid()
}

clean.by.iqr <- function(x)
{
  sx <- sort(x)
  q2 <- median(sx)
  q1 <- median(sx[sx <= q2])
  q3 <- median(sx[sx >= q2])
  iqr <- q3 - q1
  x[(q1 - 1.5 * iqr < x) & (x < q3 + 1.5 * iqr)]
}

# y <- rnorm(1000)
# qq.norm.intervals(y)