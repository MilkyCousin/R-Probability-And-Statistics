workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

alpha <- 0.05
d <- 1

X <- read.table("haker.txt", header = T)
X.0 <- data.matrix(X[colnames(X)[1]])
N <- length(X.0)

# H0: Хакерської атаки немає | Exp
# H1: Була хакерська атака | not Exp

h <- hist(X.0, probability = T, breaks = 5,
          main = "Гістограма відносних частот")

K <- length(h$breaks)
O <- h$counts
delta <- diff(h$breaks)[1]

halves <- h$breaks[-K] + delta/2
X.0.mean <- sum(O * halves)/sum(O)

l.est <- 1/X.0.mean

curve(dexp(x, l.est), col = 'red', add = T)
p.boundaries <- pexp(h$breaks, rate = l.est)
p.diff <- diff(p.boundaries[-K])
p <- c(p.diff, 1 - p.boundaries[K-1])

df.fix <- (K - 1) - d - 1

E <- N * p
chisq.stat <- sum((O - E)^2/E)
chisq.quan <- qchisq(1 - alpha, df = df.fix)
print(paste("chisq stat: ", chisq.stat))
print(paste("chisq.quan: ", chisq.quan))
print(1 - pchisq(chisq.stat, df = df.fix))