path <- "C:\\Users\\dango\\OneDrive\\Робочий стіл\\StatMagData\\CompStat"

source(paste(path, "lab1", "qqplotinterval.R", sep="\\"))

cur.data <- read.csv(paste(path, "merged.csv", sep="\\"),
                     row.names = "X")

number.of.test.sessions <- 20
# target company - cl (Colgate-Palmolive Company)
x.clo.lagged <- cur.data[-nrow(cur.data),-1]
x.clo.lagged$cl <- cur.data$cl[-1]
x.clo.lagged <- x.clo.lagged[1:(nrow(x.clo.lagged)-number.of.test.sessions),]

# first part

# full data
lm.full <- lm(cl ~ ., data=x.clo.lagged)
print("full lm")
print(summary(lm.full))

# reduced data
lm.corr <- lm(cl ~ .-clf-cmcsa-cnp, data=x.clo.lagged)
print("corr lm")
print(summary(lm.corr))

U.studentized <- rstandard(lm.corr)
hist(U.studentized, main="Histogram of studentized residues", xlab="Residue")
qq.norm.intervals(U.studentized)

plot(fitted.values(lm.corr), x.clo.lagged$cl, main="Prediction-Response diagram",
     xlab="Prediction", ylab="Response", cex=0.5)
abline(0,1,col='red')

plot(fitted.values(lm.corr), residuals(lm.corr), main="Prediction-Residue diagram",
     xlab="Prediction", ylab="Residue", cex=0.5)
abline(h=0, col='red')