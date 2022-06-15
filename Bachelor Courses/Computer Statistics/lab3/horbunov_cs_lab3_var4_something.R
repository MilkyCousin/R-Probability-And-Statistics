global.seed <- 1 #.89

rgaussmixt <- function(n, mu.1, mu.2, s.1, s.2, p)
{
  ind <- sample(c(1,2), n, replace = T, prob = c(p, 1 - p))
  u.sample <- rnorm(n, mean = c(mu.1, mu.2)[ind], sd = c(s.1, s.2)[ind])
}

norm.2.moment <- function(mu, sigma)
{
  # M[Y^2], Y ~ N(mu, sigma^2)
  mu^2 + sigma^2
}

norm.4.moment <- function(mu, sigma)
{
  # M[Y^4], Y ~ N(mu, sigma^2)
  mu^4 + 6*mu^2*sigma^2 + 3*sigma^4
}

est.mm <- function(x, mu.1, mu.2, sigma.2, p)
{
  # p*f_{N(1, theta^2)} + (1-p)*f_{N(0, 0.75^2)} => MM for theta^2
  mean(x^2)/p - (1-p)/p * (sigma.2^2 + mu.2^2) - mu.1^2
}

asympt.var <- function(theta, mu.1, mu.2, sigma.2, p)
{
  # Коефіцієнт розсіювання за теоретичною формулою
  # D[h(X)]/(H'(theta))^2
  m.x.4 <- norm.4.moment(mu.1, theta) * p + norm.4.moment(mu.2, sigma.2) * (1 - p)
  m.x.2.sq <- (p * norm.2.moment(mu.1, theta) + (1 - p) * norm.2.moment(mu.2, sigma.2))^2
  (m.x.4 - m.x.2.sq)/p^2
}

generate.estimates <- function(theta, mu.1, mu.2, sigma.2, p, n, B)
{
  set.seed(global.seed)
  # Коефіцієнт розсіювання за допомогою генерування вибірок із суміші 
  # з компонентами N(1, theta^2), N(0, 0.75^2) [тут theta - відоме]
  estimates <- numeric(B)
  for(j in 1:B)
  {
    x.mxt <- rgaussmixt(n, mu.1, mu.2, theta, sigma.2, p)
    estimates[j] <- est.mm(x.mxt, mu.1, mu.2, sigma.2, p)
  }
  estimates
}

combine.mm.var.calc <- function(theta, mu.1, mu.2, sigma.2, p, n, B, plot = T)
{
  # Обчислення коеф. розсіювання двома способами
  estimates <- generate.estimates(theta, mu.1, mu.2, sigma.2, p, n, B)
  theor.val <- asympt.var(theta, mu.1, mu.2, sigma.2, p)
  
  m.e <- mean(estimates)
  v.e <- var(estimates)
  
  obtnd.val <- n * v.e
  obtnd.bias <- sqrt(n) * (m.e - theta^2)
  
  if(plot)
  {
    # histogram of abs freq
    hist.data <- hist(estimates, plot = F)
    delta <- diff(hist.data$breaks)[1]
    
    peak <- dnorm(m.e, mean = m.e, sd = v.e^0.5) * B * delta
    
    plot(hist.data, ylim = c(0, max(peak, hist.data$counts)), 
         density = 15, main = paste("Histogram for MM-Estimates, n = ", n))
    curve(dnorm(x, mean = m.e, sd = v.e^0.5) * B * delta, 
          add = T, col = 2, lwd = 3)
    
    #hist.data <- hist(sqrt(n) * (estimates - theta^2), plot = F)
    #delta <- diff(hist.data$breaks)[1]
    
    #peak <- dnorm(m.e, mean = m.e, sd = theor.val^0.5) * B * delta
    
    #plot(hist.data, ylim = c(0, max(peak, hist.data$counts)), 
    #     density = 15, main = "Histogram for sqrt(n)*(x_n - theta^2)")
    #curve(dnorm(x, mean = 0, sd = theor.val^0.5) * B * delta, 
    #      add = T, col = 7, lwd = 3)
    
    # qq - plot for mm-estimates
    qqnorm(estimates, cex = 0.5, 
           main = paste("Normal Q-Q Plot for MM-Estimates, n = ", n))
    qqline(estimates, col = 4, lwd = 2)
  }
  
  print(paste("n:", n))
  print(paste("theoretical:", theor.val))
  print(paste("obtained:", obtnd.val))
  print(paste("bias:", obtnd.bias))
  
  list(
    est.array = estimates,
    th.v = theor.val, 
    ob.val = obtnd.val,
    ob.bias = obtnd.bias
  )
}

# Тимчасова функція, просто я морально не готовий переписувати частину коду
asympt.var.bicycle <- function(theta.sq, mu.1, mu.2, sigma.2, p)
{
  # Коефіцієнт розсіювання за теоретичною формулою
  # D[h(X)]/(H'(theta))^2
  m.x.4 <- (mu.1^4 + 6*mu.1^2*theta.sq + 3*theta.sq^2) * p + norm.4.moment(mu.2, sigma.2) * (1 - p)
  m.x.2.sq <- (p * (mu.1^2 + theta.sq) + (1 - p) * norm.2.moment(mu.2, sigma.2))^2
  (m.x.4 - m.x.2.sq)/p^2
}

# Підрахунок асимптотичного інтервалу за вибіркою та заданим alpha
conf.interval.mm <- function(x, mu.1, mu.2, s.2, p, alpha = 0.05)
{
  mm.estimate <- est.mm(x, mu.1, mu.2, s.2, p)
  quant.norm <- qnorm(1 - alpha/2)
  asympt.v <- asympt.var.bicycle(mm.estimate, mu.1, mu.2, s.2, p)
  c.interval <- mm.estimate + c(-1,1) * quant.norm * sqrt(
    asympt.v/length(x)
  )
  c.interval
}

#####

given.mu.1 <- 1
given.mu.2 <- 0
given.sigma.2 <- 0.75
given.p <- 0.6
true.theta <- 0.05

v.array <- c()
n.array <- c(100, 250, 500, 1000, 2000, 5000, 10000)
for(N in n.array)
{
  v.array <- c(v.array, combine.mm.var.calc(
    true.theta, given.mu.1, given.mu.2, given.sigma.2, given.p, 
    n = N, B = 1000, plot = T
  )$ob.val)
}
plot(n.array, v.array, type = "l", xlab = "n", ylab = "var", log = "xy")
abline(a = log(asympt.var(
  true.theta, given.mu.1, given.mu.2, given.sigma.2, given.p), 10), 
  b = 0, col = 2, lwd = 2
)
grid()

# confidence interval + bootstrap

N <- 2000
B <- 1000

counts <- c()

for(i in 1:B)
{
  u.mxt <- rgaussmixt(
    N, given.mu.1, given.mu.2, true.theta, given.sigma.2, given.p
  )
  intervals <- conf.interval.mm(
    u.mxt, given.mu.1, given.mu.2, given.sigma.2, given.p, alpha = 0.05
  )
  counts <- c(counts, 
              intervals[1] < true.theta^2 && true.theta^2 < intervals[2])
}

print(1 - mean(counts))
