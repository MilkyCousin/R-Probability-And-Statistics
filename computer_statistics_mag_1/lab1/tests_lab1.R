# Прогнозування на тестових даних. Коефіцієнти тут - оцінки МНК для 
# трьох випадків. Реалізація  брудна. 
path <- "C:\\Users\\dango\\OneDrive\\Робочий стіл\\StatMagData\\CompStat"

cur.data <- read.csv(paste(path, "merged.csv", sep="\\"),
                     row.names = "X")

number.of.test.sessions <- 20
# target company - cl (Colgate-Palmolive Company)
x.clo.lagged <- cur.data[-nrow(cur.data),-1]
x.clo.lagged$cl <- cur.data$cl[-1]
rows.to.choose <- (nrow(x.clo.lagged)-number.of.test.sessions+1):nrow(x.clo.lagged)
Y <- x.clo.lagged$cl[rows.to.choose]
X <- cbind(rep(1, number.of.test.sessions),x.clo.lagged[rows.to.choose,-1])
colnames(X)[1] <- "Intercept"


b.full <- matrix(
  c(4.43348655, 0.44885828, -0.33459201, 0.05814920, 
    -0.01041744, 0.04052351, 0.88137057), ncol=1)
# b0 + clx + cma + cme + cmg + cmi + cms

b.50.contains <- matrix(
  c(0.84364645, 0.53042112, -0.08182293), ncol=1)
# cnp + clx + cme

b.50.removed <- matrix(
  c(0.76033973, -0.07264966), ncol=1)
# clx cme

X.full <- data.matrix(X[,c("Intercept", "clx", "cma", "cme", "cmg", "cmi", "cms"),-1])
y.hat.full <- X.full%*%b.full

X.50.contains <- data.matrix(X[,c("cnp", "clx", "cme"),-1])
y.50.contains <- X.50.contains%*%b.50.contains

X.50.removed <- data.matrix(X[,c("clx", "cme"),-1])
y.50.removed <- X.50.removed%*%b.50.removed
