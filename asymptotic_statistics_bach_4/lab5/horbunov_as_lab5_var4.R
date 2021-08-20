library(corrplot)
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

Y <- wine.data$Site[wine.data$Site != 2] 

Y[Y == 3] <- 0

wine.data <- wine.data[
  wine.data$Site != 2,
  c("Proanthocyanins", "Color_intensity", "Hue", "OD")
  ]

cor.p <- cor(wine.data)
cor.s <- cor(wine.data, method = "spearman")

par(mfrow = c(1,2))
corrplot(cor.p, method = "color", addCoef.col = "black", 
         main = "Pearson Correlation", mar=c(5,0,8,1))
corrplot(cor.s, method = "color", addCoef.col = "black", 
         main = "Spearman Correlation", mar=c(5,0,8,1))
par(mfrow=c(1,1))

alpha <- 0.05

logit.regr.1 <- glm(Y ~ Proanthocyanins + Hue,
                  data = wine.data, family = "binomial")
print(summary(logit.regr.1))
Y.pred.1 <- ifelse(logit.regr.1$fitted.values > 0.5, 1, 0)

print(
  confint(logit.regr.1, level = 1 - alpha)
)

print(
  classification_summary.binom(Y, Y.pred.1)
)

stop()

logit.regr.0 <- glm(Y ~ Proline + Malic_acid + Alcogol,
                    data = wine.data, family = "binomial")
print(summary(logit.regr.0))
Y.pred.0 <- ifelse(logit.regr.0$fitted.values > 0.5, 1, 0)

print(
  confint(logit.regr.0, level = 1 - alpha)
)

print(
  classification_summary.binom(Y, Y.pred.0)
)


