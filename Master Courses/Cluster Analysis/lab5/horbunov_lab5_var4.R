library('rgl')
options(rgl.printRglwidget = TRUE)

data <- read.table("v4.txt", header=T)
plot(data)

mfrow3d(2, 2)
plot3d(data[,c(2,1,5)], main="y,x,v")
plot3d(data[,c(2,3,5)], main="y,z,v")
plot3d(data[,c(2,4,5)], main="y,u,v")
plot3d(data[,c(2,4,5)], main="y,u,v")
mfrow3d(1, 1)

library('mclust')

# Model Fit
data.mclust <- Mclust(data=data)
plot(data.mclust)

pal <- c('red', 'green', 'blue', 'orange')
plot3d(data[,c(2,4,5)], col=pal[data.mclust$classification])

# Dimension Reduction
data.mclust.dr <- MclustDR(object=data.mclust, lambda=0.5)
plot(data.mclust.dr)

dr.b <- data.mclust.dr$basis
plot3d((data.matrix(data)%*%dr.b)[,1:3], col=pal[data.mclust$classification])

par(mfrow=c(2,2))
# Dimensionality Per Cluster
cov.estimates <- data.mclust$parameters$variance$sigma
for(k in 1:4)
{
  pca.k <- princomp(covmat=cov.estimates[,,k])
  plot(pca.k, main=paste("Cluster #", k, sep=""))
  print(summary(pca.k))
}
par(mfrow=c(1,1))

# Other

l.pca.1 <- matrix(loadings(princomp(covmat=cov.estimates[,,1]))[,1], ncol=1)
l.pca.2 <- loadings(princomp(covmat=cov.estimates[,,2]))[,1:2]
l.pca.3 <- loadings(princomp(covmat=cov.estimates[,,3]))[,1:2]
l.pca.4 <- loadings(princomp(covmat=cov.estimates[,,4]))[,1:2]

print(t(l.pca.1)%*%l.pca.2[,1])
print(t(l.pca.1)%*%l.pca.2[,2])

print(t(l.pca.1)%*%l.pca.3[,1])
print(t(l.pca.1)%*%l.pca.3[,2])

print(t(l.pca.1)%*%l.pca.4[,1])
print(t(l.pca.1)%*%l.pca.4[,2])