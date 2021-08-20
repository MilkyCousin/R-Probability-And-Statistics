workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()


distr.seq <- read.table("distr4.txt", header = TRUE)$X

l <- length(distr.seq)
s <- sqrt(var(distr.seq))
m <- mean(distr.seq)

brk <- 40 #seq(min(distr.seq), max(distr.seq))  #((0:48))*0.1

plot.hist <- function(t, distrib.expr, plt.title = "")
{
  hist(t, probability = TRUE, breaks = brk, main = plt.title, ylim=c(0,0.2))
  curve(distrib.expr(x), add=T, lwd=2, col='red')
}

# norm
plot.hist(distr.seq, function(x) {dnorm(x, mean = m, sd = s)}, 
          "Histogram + Estimated Normal Distribution")

# exp(norm)
lnorm.var <- log((s/m)^2 + 1)
lnorm.mean <- log(m/sqrt(lnorm.var))
plot.hist(distr.seq, function(x) {dlnorm(x, meanlog = lnorm.mean, sdlog = sqrt(lnorm.var))}, 
          "Histogram + Estimated Lognormal Distribution")

# exp
exp.lambda <- (l-1)/l * 1/m
plot.hist(distr.seq, function(x) {dexp(x, rate = exp.lambda)}, 
          "Histogram + Estimated Exponential Distribution")

# chisq
plot.hist(distr.seq, function(x) {dchisq(x, df=m)},
          "Histogram + Estimated Chi-Square Distribution")

# QQ,PP-diagrams
ecdf.vals <- ecdf(distr.seq)

plot.pp <- function(
  t, p.distr.expr,
  x.lab = "Theoretical P", y.lab = "Empirical P", plt.title = ""
)
{
  l <- length(t)
  st <- sort(t)
  plot(p.distr.expr(st), (1:l)/l, cex=0.25, asp = 1,
       ylab=y.lab, xlab=x.lab, main=plt.title)
  abline(0,1,col=2)
}

plot.pp(distr.seq, function(x) {pexp(sort(x), rate = exp.lambda)},
        plt.title = "Exp-ECDF comparison")

plot.pp(distr.seq, function(x) {pchisq(sort(x), df=m)},
        plt.title = "ChiSq-ECDF comparison")

plot(qexp(ppoints(distr.seq), rate=exp.lambda), sort(distr.seq), cex = 0.25, 
     xlim = c(0, 6), ylim = c(0, 6), main = "Exp-ECDF comparison",
     xlab = "Theoretical Q", ylab = "Empirical Q")
abline(lm(
  sort(distr.seq)~qexp(ppoints(distr.seq), rate=exp.lambda)
  )$coefficients)

# QQ with confidence intervals
qqexp.conf <- function(x, rate.value = 1, alpha = .05, k = 10^3)
{
  l <- length(x)
  sx <- sort(x)
  
  generated <- matrix(
    rexp((k * l), rate = rate.value), 
    nrow = l, ncol = k
    )
  
  generated <- apply(generated, 2, sort)
  
  upper <- apply(generated, 1, quantile, probs = 1 - alpha * .5)
  lower <- apply(generated, 1, quantile, probs = alpha * .5)
  
  exp.quant <- qexp((1:l-.5)/l, rate = rate.value)
  
  plot(
    c(exp.quant, exp.quant, exp.quant),
    c(upper, lower, sx), type = "n",
    xlab = "Theoretical quantiles",
    ylab = "Empirical quantiles",
    main = paste("Q-Q plot with confidence intervals, alpha = ", alpha)
  )
  
  points(exp.quant, sx, col = "orange", cex=0.45)
  segments(exp.quant, lower, exp.quant, upper, col = "purple")
  
  abline(0, 1, col = "black", lwd=1)
}

qqexp.conf(distr.seq, 1/m)