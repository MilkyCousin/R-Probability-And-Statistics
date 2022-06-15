workspace_setup <- function(seed.number = 0)
{
  set.seed(seed.number)
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

distr.seq <- read.table("distr4.txt", header = TRUE)$X

clean.by.iqr <- function(x)
{
  sx <- sort(x)
  q2 <- median(sx)
  q1 <- median(sx[sx <= q2])
  q3 <- median(sx[sx >= q2])
  iqr <- q3 - q1
  x[(q1 - 1.5 * iqr < x) & (x < q3 + 1.5 * iqr)]
}

l <- length(distr.seq)
s <- sqrt(var(distr.seq))
m <- mean(distr.seq)

exp.lambda <- (l-1)/l * 1/m

histogram.confidence <- function(
  x, r = 1, alpha = .05, k = 10^3, eps = .5, step = .25, y.lim=c(0,1)
  ) # suppose we have exp distribution
{
  n <- length(x)
  exp.generated <- matrix(rexp((n * k), rate = r), nrow = n, ncol = k)
  
  top.margin <- max(apply(exp.generated, 2, 
                          function(t) {max(clean.by.iqr(t))}), max(x))
  low.margin <- min(apply(exp.generated, 2, 
                          function(t) {min(t)}), min(x))
  
  brk.margins <- seq(low.margin - eps, top.margin + eps, step)

  hist.data <- apply(
    exp.generated, 2, 
    function(t) {
      h <- hist(clean.by.iqr(t), breaks = brk.margins, plot=F); 
       (
        function(u) {v <- numeric(length(brk.margins)); v[1:length(u)] = u; v}
       )((h$counts/diff(h$breaks))/sum(h$counts))
      }
    )
  
  q.low <- apply(hist.data, 1, quantile, probs = alpha * .5)
  q.top <- apply(hist.data, 1, quantile, probs = 1 - alpha * .5)
  
  hst <- hist(x, breaks = brk.margins, probability = T, ylim = y.lim,
              main = paste(
                "Histogram with intervals for relative frequencies, alpha = ", 
                alpha
                ))
  curve(dexp(x, rate = r), add = T, lwd = 2, col = 'red')
  segments(brk.margins + step * .5, q.low, 
           brk.margins + step * .5, q.top, 
           lwd=2, col='blue')
  legend("topright", 
         legend = c(
           "Probability density function curve", "Confidence intervals"),
         col = c("red", "blue"), lwd = c(2,2)
         )
}