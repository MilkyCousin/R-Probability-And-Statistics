library(MASS)

CODENAMES <-c("S. glaber, Norway", 
              "S. linnarssoni, Norway", 
              "S. linnarssoni, Sweden")
RGB <- c("Red", "Green", "Blue")
GO <- c("Green", "Orange")

tri1 <- read.table("glaber.dat", header = T, row.names = ".")
tri2 <- read.table("linnars1.dat", header = T, row.names = ".")
tri3 <- read.table("linnars2.dat", header = T, row.names = ".")

n1 <- nrow(tri1)
n2 <- nrow(tri2)
n3 <- nrow(tri3)

tri <- cbind(
  rbind(tri1, tri2, tri3),
  rep(1:3, c(n1, n2, n3)), 
  c(rep("G", n1), rep("L", n2 + n3))
)
colnames(tri)[c(5,6)] <- c("K", "Code")

tri.manova <- manova(cbind(tri$W1, tri$W2, tri$W3, tri$L2) ~ tri$Code)

tri.s <- tri[tri$Code == "L",]
tri.manova.local <- manova(cbind(tri.s$W1, tri.s$W2, tri.s$W3, tri.s$L2) ~ tri.s$K)

par(mfrow = c(1,4))
boxplot(W1 ~ Code, data = tri)
boxplot(W2 ~ Code, data = tri)
boxplot(W3 ~ Code, data = tri)
boxplot(L2 ~ Code, data = tri)
par(mfrow = c(1,1))

plot(tri[,1:4], col = RGB[tri$K], cex= 0.75)
plot(tri[,1:4], col = GO[(tri$Code == "L") + 1], cex= 0.75)

# Підрахунок статистики тесту Чоу
lm.h0 <- lm(L2  ~ W1, data = tri)
resid.h0 <- residuals(lm.h0)

lm.h1.subset1 <- lm(L2  ~ W1, data = subset(tri, Code == "G"))
lm.h1.subset2 <- lm(L2  ~ W1, data = subset(tri, Code == "L"))
resid.h1 <- c(residuals(lm.h1.subset1), residuals(lm.h1.subset2))

N <- nrow(tri)
alpha <- 0.05
# Поріг тесту
F.theor <- qf(1 - alpha, 2, N - 4)
Sh0 <- sum(resid.h0^2)
Sh1 <- sum(resid.h1^2)
# Статистика тесту
F.emp <- ((1 / 2) * (Sh0 - Sh1)) / ((1 / (N - 4)) * Sh1)
#

tri.pca <- princomp(tri[,1:4], cor = F)
summary(tri.pca)

plot(tri.pca, main = "Діаграма власних чисел")

# Виводимо матрицю навантажень
tri.pca$loadings

barplot(t(tri.pca$loadings[,1]), main = "Навантаження першої компоненти")
barplot(t(tri.pca$loadings[,2]), main = "Навантаження другої компоненти")
barplot(t(tri.pca$loadings[,3]), main = "Навантаження третьої компоненти")
barplot(t(tri.pca$loadings[,4]), main = "Навантаження четвертої компоненти")

plot(tri.pca$scores[,1], tri.pca$scores[,2],
     xlab = "Comp.1", ylab = "Comp.2", col = RGB[tri$K], cex= 0.75, lwd = 2)
grid()
legend("topright", legend = CODENAMES, col = RGB, pch = 1)

tri.proj <- data.frame(cbind(tri.pca$scores[,c(1,2)], tri$Code))
colnames(tri.proj)[3] <- "Code"

tri.pca.lda <- lda(Code ~ ., data = tri.proj)
tri.pca.lda

tri.pca.lda.cv <- lda(Code ~ ., data = tri.proj, CV = T)
table(forecast = tri.pca.lda.cv$class, real = tri$Code)

plot(tri.pca$scores[,1], tri.pca$scores[,2],
     xlab = "Comp.1", ylab = "Comp.2", 
     col = GO[(tri$Code == "L")+1], cex= 0.75, lwd = 2)
grid()
legend("topright", legend = c("S. glaber", "S. linnarssoni"), col = GO, pch = 1)

scaling.tri <- tri.pca.lda[4]$scaling
abline(scaling.tri[c(2,1)], col = "purple", lty = 2)

legend("topright", legend = c("S. glaber", "S. linnarssoni"), col = GO, pch = 1)


ty.lda <- function(x, groups){
  x.lda <- lda(groups ~ ., as.data.frame(x))
  
  gr <- length(unique(groups))   ## groups might be factors or numeric
  v <- ncol(x) ## variables
  m <- x.lda$means ## group means

  w <- array(NA, dim = c(v, v, gr))
  
  for(i in 1:gr){
    tmp <- scale(subset(x, groups == unique(groups)[i]), scale = FALSE)
    w[,,i] <- t(tmp) %*% tmp
  }
  
  W <- w[,,1]
  for(i in 2:gr)
    W <- W + w[,,i]
  
  V <- W/(nrow(x) - gr)
  iV <- solve(V)
  
  class.funs <- matrix(NA, nrow = v + 1, ncol = gr)
  colnames(class.funs) <- paste("group", 1:gr, sep=".")
  rownames(class.funs) <- c("constant", paste("var", 1:v, sep = "."))
  
  for(i in 1:gr) {
    class.funs[1, i] <- -0.5 * t(m[i,]) %*% iV %*% (m[i,])
    class.funs[2:(v+1) ,i] <- iV %*% (m[i,])
  }
  
  x.lda$class.funs <- class.funs
  
  return(x.lda)
}


ty.lda.fit <- ty.lda(tri.proj[,c(1,2)], tri$Code)
ty.lda.fit