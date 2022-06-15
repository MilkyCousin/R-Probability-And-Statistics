path <- "C:\\Users\\dango\\OneDrive\\Робочий стіл\\StatMagData\\CompStat"

source(paste(path, "lab1", "qqplotinterval.R", sep="\\"))

cur.data <- read.csv(paste(path, "merged.csv", sep="\\"),
                     row.names = "X")

# target company - cl (Colgate-Palmolive Company)

number.of.test.sessions <- 20
x.clo.lagged <- cur.data[-nrow(cur.data),-1]
x.clo.lagged$cl <- cur.data$cl[-1]
x.clo.lagged <- x.clo.lagged[1:(nrow(x.clo.lagged)-number.of.test.sessions),]

n.sessions <- 50
x.nrow <- nrow(x.clo.lagged)

x.recent <- x.clo.lagged[(x.nrow-n.sessions+1):x.nrow,]
x.nrow.rec <- nrow(x.recent)

lm.partial <- lm(cl ~ ., data = x.recent)
print(summary(lm.partial))

remained_cols <- c("cl", "cnp", "clx", "cme")
x.remained <- x.recent[,remained_cols]
rownames(x.remained) <- 1:nrow(x.remained)

lm.remained <- lm(cl ~ .-1, data = x.remained)
print(summary(lm.remained))

# analysis of residuals

U.studentized <- rstandard(lm.remained)

hist(U.studentized, prob = T,
     main = "Histogram of studentized residuals",
     xlab = "Studentized residue", ylab = "Density")
curve(dnorm(x), col='red', add=T)
grid()

qq.norm.intervals(U.studentized)

# finally, influence plot

plot(fitted.values(lm.remained), x.remained$cl, main="Prediction-Response diagram",
     xlab="Prediction", ylab="Response", cex=0.5)
abline(0,1,col='red')

plot(fitted.values(lm.remained), residuals(lm.remained), main="Prediction-Residue diagram",
     xlab="Prediction", ylab="Residue", cex=0.5)
abline(h=0, col='red')

lm.influence <- car::influencePlot(lm.remained)
grid()

# after influence diagnostics
x.removed <- x.remained[-c(1,19,34,35),]
lm.removed <- lm(cl ~ .-cnp-1, data = x.removed)
print(summary(lm.removed))