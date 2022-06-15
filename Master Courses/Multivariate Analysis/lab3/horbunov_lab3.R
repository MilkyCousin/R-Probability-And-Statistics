library(readxl)
library(pvclust)
library(rgl)
options(rgl.printRglwidget = TRUE)

birds.dat <- read_xls("./birds.xls")
birds.dat <- data.frame(birds.dat)[,-1]

attach(birds.dat)

# Нормування ваги
W <- w^(1/3)

# Зводимо до єдиного виміру (з дюймів до міліметрів)
TL <- log(tl)
AE <- log(ae)
LBH <- log(lbh)
LH <- log(lh^25.4)
LF <- log(lf^25.4)
LT <- log(lt^25.4)
WS <- log(ws^25.4)
LKS <- log(lks^25.4)

# Побудова таблиці
data1 <- cbind(W, TL, AE, LBH, LH, LF, LT, WS, LKS)

pv1 <- pvclust(data1, method.hclust = "average", method.dist = "abscor")
plot(pv1)
pvrect(pv1, lty = 2)

clusters.picked <- pvpick(pv1)
M <- length(clusters.picked$clusters)

color.choose <- function(features, pos0 = 1)
{
  sapply(features, function(feature) {
    for(m in 1:M)
    {
      if(feature%in%clusters.picked$clusters[[m]])
      {
        return(pos0 + m)
      }
    }
    return(pos0)
  })
}

FEATURES <- colnames(data1)
colored <- color.choose(FEATURES)

# MDS 2d
mds2d <- cmdscale(as.dist(1 - abs(cor(data1))), k = 2)
plot(mds2d[,1], mds2d[,2], type = "n", xlab = "Axis1", ylab = "Axis2",
     main = "Multi-dimensional scaling, 2d")
grid()
text(mds2d[,1], mds2d[,2], labels = FEATURES)#, col = colored)

# MDS 3d
mds3d <- cmdscale(as.dist(1 - abs(cor(data1))), k = 3)
plot3d(mds3d[,1], mds3d[,2], mds3d[,3], type = "n", 
       xlab = "Axis1", ylab = "Axis2", zlab = "Axis3",
       main = "Multi-dimensional scaling, 3d")
text3d(mds3d[,1], mds3d[,2], mds3d[,3], texts = FEATURES)#, col = colored)

# PCA
pca <- princomp(data1, cor = T)
summary(pca)
plot(pca)
plot(pca$loadings[,1], pca$loadings[,2], type = "n",
     xlab = "Comp.1", ylab = "Comp.2", main = "PC loadings")
grid()
text(pca$loadings[,1], pca$loadings[,2], label = FEATURES)#, col = colored)

plot3d(pca$loadings[,1], pca$loadings[,2], pca$loadings[,3], type = "n", 
       xlab = "Comp.1", ylab = "Comp.2", zlab = "Comp.3",
       main = "PC loadings")
text3d(pca$loadings[,1], pca$loadings[,2], pca$loadings[,3], 
       texts = FEATURES)#, col = colored)