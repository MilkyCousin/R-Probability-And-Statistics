set.seed(0)

n <- 300
x <- rnorm(n, 0, 1)
e <- rnorm(n, 0, sqrt(0.5))
g <- function(t) { 2 * abs(t) }
y <- g(x) + e

plot(y ~ x); grid()

moving.average <- function(x, y, h)
{
  ma.univar <- function(t) 
  { 
    ones <- (abs(t - x) < h / 2)
    sum(y * ones) / sum(ones) 
  }
  function(t) sapply(t, ma.univar)
}

moving.median <- function(x, y, h)
{
  mm.univar <- function(t) 
  { 
    median(y[abs(t - x) < h / 2])
  }
  function(t) sapply(t, mm.univar)
}

loc.lin.regr <- function(x, y, h, K, corr.w = 1)
{
  llr.univar <- function(x0)
  {
    w <- K((x0 - x) / h)
    loc.lm <- lm(y ~ I(x - x0), weights = corr.w * w)
    coef(loc.lm)[1]
  }
  function(t) sapply(t, llr.univar)
}

K.e <- function(t) { 0.75 * (1 - t^2) * (abs(t) < 1) }
K.b <- function(t) { (1 - t^2)^2 * (abs(t) < 1) }

#x <- c(x, rep(c(-2, 0, 2), 5))
#y <- c(y, rep(c(10, -5, 10), 5))

plot(y ~ x); grid()


H <- seq(0.01, 1, 0.01)
# ma - h = 0.75

# Moving average

for(j in 1:length(H))
{
  png(file = paste(
    "/home/coolpenguin/R/nonparam/lab5/ma/ma", formatC(j, width = floor(log(length(H), base=10)) + 1, flag = "0"),
    ".png", sep = ""
    ))
  plot(y ~ x, main = paste("Moving Average, h =", H[j])); grid()
  curve((moving.average(x, y, H[j]))(t), xname = "t", col = 2, add = T, type = "s")
  dev.off()
}

# Moving median

for(j in 1:length(H))
{
  png(file = paste(
    "/home/coolpenguin/R/nonparam/lab5/mm/mm", formatC(j, width = floor(log(length(H), base=10)) + 1, flag = "0"),
    ".png", sep = ""
  ))
  plot(y ~ x, main = paste("Moving Median, h =", H[j])); grid()
  curve((moving.median(x, y, H[j]))(t), xname = "t", col = 2, add = T, type = "s")
  dev.off()
}

# loc-lin regr, epan kern
for(j in 1:length(H))
{
  png(file = paste(
    "/home/coolpenguin/R/nonparam/lab5/ll/ll", formatC(j, width = floor(log(length(H), base=10)) + 1, flag = "0"),
    ".png", sep = ""
  ))
  plot(y ~ x, main = paste("Loc-Linear Regression, h =", H[j])); grid()
  curve((loc.lin.regr(x, y, H[j], K.e))(t), xname = "t", col = 2, add = T, type = "s")
  dev.off()
}

plot(y ~ x); grid()
curve((moving.average(x, y, 0.69))(t), xname = "t", col = 2, add = T, type = "s")
curve((moving.median(x, y, 0.69))(t), xname = "t", col = 4, add = T, type = "s")
curve((loc.lin.regr(x, y, 0.6, K.e))(t), xname = "t", col = 6, add = T, type = "s")

legend("bottomleft", legend = c("Moving Average", "Moving Median", "Loc-Linear Regression"), lwd = 1, col = 2*(1:3))

# corrected

U <- y - (loc.lin.regr(x, y, 0.4, K.e))(x)
mu <- median(abs(U))
deltaj <- K.b(U / (6 * mu))

plot(y ~ x); grid()

curve((loc.lin.regr(x, y, 0.4, K.e))(t), xname = "t", col = 1, add = T, type = "s")
curve((loc.lin.regr(x, y, 0.4, K.e, corr.w = deltaj))(t), xname = "t", col = 2, add = T, type = "s")

for(j in c(3,4))
{
  
  U <- y - (loc.lin.regr(x, y, 0.4, K.e, corr.w = deltaj))(x)
  mu <- median(abs(U))
  deltaj <- K.b(U / (6 * mu))
  curve((loc.lin.regr(x, y, 0.4, K.e, corr.w = deltaj))(t), xname = "t", col = j, add = T, type = "s")
}

legend("bottomleft", 
       legend = c("До виправлення", "Перше виправлення", "Друге виправлення", "Третє виправлення"), lwd = 1, col = 1:4)