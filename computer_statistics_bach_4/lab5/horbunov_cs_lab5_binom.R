# H0 : Binom(0.5, 2), p := 0.5
# H1 : Binom(0.6, 2), q := 0.6
# n = 65

set.seed(0)

log.lr.test <- function(x, m, p, q, alpha = 0.05, to.print = F)
{
  x.sum <- sum(x)
  n <- length(x)
  c.alpha <- qbinom(1 - alpha, n * m, p)
  h.res <- c.alpha >= x.sum
  if(to.print)
  {
    print(paste("sum(X) = ", x.sum, 
                ifelse(h.res, "<=", ">"),
                c.alpha, "= c_{alpha}")) 
  }
  list(hypothesis = 1 - h.res, statistic = x.sum) # 0 - H0, 1 - H1
}

p <- 0.5
q <- 0.6
m <- 2
n.fixed <- 65

# Приклад застосування
print(log.lr.test(rbinom(65, 2, 0.5), 2, 0.5, 0.6, to.print = T))

# Імітаційне моделювання
B <- 10000

# Генеруємо вибірки за розподілом при виконанні нульової гіпотези
# У відповідні масиви вносимо результат тесту та значення статистики
counts.0.1 <- numeric(B)
x.stat.0 <- numeric(B)
for(j in 1:B)
{
  x.0 <- rbinom(n.fixed, m, p)
  r <- log.lr.test(x.0, m, p, q)
  counts.0.1[j] <- r$hypothesis
  x.stat.0[j] <- r$statistic
}
# Оцінка для ймовірності помилки першого роду
print(
  paste("Estimated probability of I type error:", 
        mean(counts.0.1))
  )

# Генеруємо вибірки за розподілом при виконанні альтернативної гіпотези
# У відповідні масиви вносимо результат тесту та значення статистики
counts.1.1 <- numeric(B)
x.stat.1 <- numeric(B)
for(j in 1:B)
{
  x.1 <- rbinom(n.fixed, m, q)
  r <- log.lr.test(x.1, m, p, q)
  counts.1.1[j] <- r$hypothesis
  x.stat.1[j] <- r$statistic
}
# Оцінка для ймовірності помилки другого роду
print(
  paste("Estimated probability of II type error:", 
        1 - mean(counts.1.1))
  )

# Гістограма статистик в залежності від параметра ймовірності
min.stat <- min(x.stat.0, x.stat.1)
max.stat <- max(x.stat.0, x.stat.1)

hist(x.stat.0, col = 'red', xlim = c(min.stat, max.stat), ylim = c(0, 0.08),
     freq = F, breaks = 20, main = "Histogram of LR-statistics", angle = -45,
     density = 15, xlab = "Value", panel.first = grid())
hist(x.stat.1, col = 'blue', xlim = c(min.stat, max.stat), ylim = c(0, 0.08),
     freq = F, breaks = 20, add = T, density = 15, angle = 45)

# Емпірична оцінка порогу тесту
c.b <- quantile(x.stat.0, 1 - 0.05)
print(c.b)
abline(v = c.b, lwd = 2, col = 'darkgreen')

# Знаходження оптимального обсягу вибірки
