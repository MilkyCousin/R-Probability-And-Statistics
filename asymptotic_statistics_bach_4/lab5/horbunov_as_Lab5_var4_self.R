set.seed(0)

workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

classification_summary.binom <- function(y.real, y.pred)
{
  confusion.matrix <- table(y.pred, y.real)
  
  PRE <- confusion.matrix[1,1]/(confusion.matrix[1,1] + confusion.matrix[1,2])
  REC <- confusion.matrix[1,1]/(confusion.matrix[1,1] + confusion.matrix[2,1])
  
  ACC <- sum(diag(confusion.matrix))/sum(confusion.matrix)
  
  list(
    conf.matr = confusion.matrix,
    precision = PRE,
    recall    = REC,
    accuracy  = ACC
  )
}

g <- function(t)
{
  1/(1 + exp(-t))
}

f <- function(b.0, b, x)
{
  g(b.0 + t(b)%*%x)
}

lik.function <- function(b.0, b, X, Y)
{
  f.0 <- function(u) {g(b.0 + t(b)%*%u)}
  f.val <- apply(X, 1, f.0)
  lik.val <- prod(f.val^Y * (1 - f.val)^(1-Y))
  log(lik.val)
}

mle.estimation <- function(X, Y, z.0 = rep(0.01, ncol(X) + 1))
{
  to.minimize <- function(z) {-lik.function(z[1], z[-1], X, Y)}
  nlm.struct <- nlm(to.minimize, z.0)
  #print(summary(nlm.struct))
  nlm.struct$estimate
}

# pract

wine.data <- read.csv2("wine.csv")

Y <- wine.data$Site[wine.data$Site%in%c(1,2)] - 1 # 0 - 1st, 1 - 2nd

wine.data <- wine.data[
  wine.data$Site%in%c(1,2),
  c("Proline", "Alcogol", "Malic_acid", "Ash")
]

res <- mle.estimation(wine.data, Y)

aposterior.proba.estim <- function(x)
{
  f(res[1], res[-1], x)
}

predictions <- apply(wine.data, 1, function(u) {aposterior.proba.estim(u)})
predicted.values <- predictions > 0.5

print(classification_summary.binom(Y, predicted.values))