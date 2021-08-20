set.seed(0)

workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

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

workspace_setup()

wine.data <- read.csv2("wine.csv")

# var 4 -> vinohradnyki #1,2; Proline; Alcohol; Malic Acid; Ash

Y <- wine.data$Site[wine.data$Site%in%c(1,2)] - 1 # 0 - 1st, 1 - 2nd

wine.data <- wine.data[
  wine.data$Site%in%c(1,2),
  c("Proline", "Alcogol", "Malic_acid", "Ash")
]

n <- length(Y)

alpha <- 0.5
n.train <- round(alpha * n)

train.idx <- sample(1:n, n.train)
X.train <- wine.data[train.idx,]
X.test <- wine.data[-train.idx,]

logit.regr.0 <- glm(Y[train.idx] ~ Proline + Malic_acid + Alcogol,
                    data = X.train, family = "binomial")

print(summary(logit.regr.0))

Y.test.pred.0 <- ifelse(predict(logit.regr.0, X.test) > 0.5, 1, 0)
print(classification_summary.binom(Y[-train.idx], Y.test.pred.0))

