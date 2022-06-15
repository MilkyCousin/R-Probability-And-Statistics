global.seed <- 2

rgaussmixt <- function(n, mu.1, mu.2, s.1, s.2, p)
{
  ind <- sample(c(1,2), n, replace = T, prob = c(p, 1 - p))
  u.sample <- rnorm(n, mean = c(mu.1, mu.2)[ind], sd = c(s.1, s.2)[ind])
}

est.mm <- function(x, m.1, m.2, s.2, p)
{
  est.r <- mean(x^2)/p - (1-p)/p * (s.2^2 + m.2^2) - m.1^2
  sqrt(abs(est.r))
}

d.gaussmixt <- function(s.1, s.2, m.1, m.2, p, x)
{
  p*dnorm(x, m.1, s.1) + (1-p)*dnorm(x, m.2, s.2)
}

ll.mxt <- function(s.1, s.2, m.1, m.2, p, x)
{
  sum(log(d.gaussmixt(s.1, s.2, m.1, m.2, p, x)))
  #prod(d.gaussmixt(s.1, s.2, m.1, m.2, p, x))
}

est.ml <- function(x, m.1, m.2, s.2, p, x.0 = est.mm(x, m.1, m.2, s.2, p))
{
  calc <- nlm(
    function(t){
      -ll.mxt(t, s.2, m.1, m.2, p, x)
    }, x.0 #, print.level = 2
    )
  t <- seq(0, 0.25, 0.001)
  m <- length(t)
  y <- numeric(m)
  for(i in 1:m)
  {
    y[i] <- ll.mxt(t[i], s.2, m.1, m.2, p, x)
  }
  plot(t, y, type='l', main='log-likelihood plot'); grid()
  abline(v = calc$estimate, col = 'red')
  calc$estimate^2
}

####

h <- function(t, m.1, m.2, s.1, s.2, p)
{
  g <- p/sqrt(2*pi) * exp(-0.5 * (t - m.1)^2/s.1^2) * ((t - m.1)^2 * s.1^(-4) - s.1^(-2))
  l <- d.gaussmixt(s.1, s.2, m.1, m.2, p, t)
  g^2/l
}

fisher.info <- function(m.1, m.2, s.1, s.2, p, a = -15, b = -a)
{
  integration.result <- integrate(
    function(x)
    {
      h(x, m.1, m.2, s.1, s.2, p)
    },
    a, b
  )
  integration.result
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
    estimates[j] <- est.ml(x.mxt, mu.1, mu.2, sigma.2, p)
  }
  estimates
}

mle.asympt.var <- function(theta, mu.1, mu.2, sigma.2, p)
{
  I <- fisher.info(mu.1, mu.2, theta, sigma.2, p)
  1/I$value
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

mm.asympt.var <- function(theta, mu.1, mu.2, sigma.2, p)
{
  # Коефіцієнт розсіювання за теоретичною формулою (для моментної оцінки)
  # D[h(X)]/(H'(theta))^2
  m.x.4 <- norm.4.moment(mu.1, theta) * p + norm.4.moment(mu.2, sigma.2) * (1 - p)
  m.x.2.sq <- (p * norm.2.moment(mu.1, theta) + (1 - p) * norm.2.moment(mu.2, sigma.2))^2
  (m.x.4 - m.x.2.sq)/p^2
}

are.mle.mm <- function(theta, mu.1, mu.2, sigma.2, p)
{
  mle.v <- mle.asympt.var(theta, mu.1, mu.2, sigma.2, p)
  gmm.v <- mm.asympt.var(theta, mu.1, mu.2, sigma.2, p)
  print(c(mle.v, gmm.v))
  mle.v/gmm.v
}

####

given.mu.1 <- 1
given.mu.2 <- 0
given.sigma.2 <- 0.75
given.p <- 0.6
true.theta <- 0.05

N <- 1000
u.mxt <- rgaussmixt(
  N, 
  given.mu.1, given.mu.2, 
  true.theta, given.sigma.2, 
  given.p
)

options(warn = -1)
estimated.ml <- est.ml(
  u.mxt, given.mu.1, given.mu.2, given.sigma.2, given.p
)
print(estimated.ml)
options(warn = 0)

###

#I.d <- fisher.info(given.mu.1, given.mu.2, true.theta, given.sigma.2, given.p)
#I <- I.d$value
#V <- 1/I

#UU <- generate.estimates(0.05, 1, 0, 0.75, 0.6, 1000, 1000)
#hist(sqrt(1000)*(UU[UU<.6] - true.theta^2), main = "Histogram of sqrt(n)*(theta_n - theta^2)", probability = T); grid()
#curve(dnorm(x, mean = sqrt(1000)*(mean(UU[UU < .6]) - true.theta^2), sd = (sqrt(1000*var(UU[UU < .6]))), col = 'red', add = T))

####

print('are')
print(are.mle.mm(true.theta, given.mu.1, given.mu.2, given.sigma.2, given.p))