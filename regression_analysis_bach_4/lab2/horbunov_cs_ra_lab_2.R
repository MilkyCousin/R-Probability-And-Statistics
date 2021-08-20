workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

df <- read.csv("ms3.csv", sep = ' ')[c('MKR1', 'MKR2')]
plot(df)
df <- df[-10,]
row.names(df) <- NULL
n <- nrow(df)
d <- 2

x <- cbind(1 + numeric(n), data.matrix(df[1]))
y <- data.matrix(df[2])

hist(x[,2], probability = T, main = "Гістограма результатів МКР1", breaks = 10)
hist(y, probability = T, main = "Гістограма результатів МКР2")
curve(dnorm(x, mean(y), sd(y)), col = 'blue', add = T)

a <- t(x)%*%x
a.inv <- solve(a)
b <- a.inv%*%t(x)%*%y

plot(df)
abline(b, col = 'red')

y.h <- x%*%b # predictions
u <- y - y.h # residuals

# descriptive statistics
hist(u, probability = T, breaks = 10, main = "Гістограма залишків")
curve(dnorm(x, mean(u), sd(u)), col = 'blue', add = T)

qqnorm(u, main = "QQ-діаграма залишків")
qqline(u)

plot(y.h, u, main = "Діаграма 'Прогноз-залишки' для початкових даних.")
plot(y.h, y, main = "Діаграма 'Прогноз-відгук' для початкових даних.")
abline(lm(y ~ y.h), col = 'red')

# r-squared
r.sq <- sum(u^2)/sum((y - mean(y))^2)
print(1 - r.sq)

# variance of errors estimation
err.var.estim <- sum(u^2)/(n - d) 
print(err.var.estim)

# variance of coefficients estimation
coef.var.estim <- err.var.estim * diag(a.inv)
print(sqrt(coef.var.estim))

# t-test
t.coef.0 <- -b/sqrt(coef.var.estim)
print(t.coef.0)

alpha <- 0.05
q.t <- qt(1 - alpha/2, df = n - d)
print(2*pt(t.coef.0, df = n - d))

err.t <- q.t * sqrt(coef.var.estim)

# c.i. for b
l.b <- b - err.t
u.b <- b + err.t
print("c.i.")
print(l.b)
print(u.b)

# F-quantile
m <- d - 1
q.f <- qf(1 - alpha, m, n - m - 1)
print(q.f)
f.emp <- (n - m - 1)/m * r.sq/(1 - r.sq)
print(f.emp)