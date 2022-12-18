library('kernlab')
library('rgl')

options(rgl.printRglwidget = TRUE)

data <- read.table("F4p.txt", header=F)
plot(data, cex=0.25)

print(diag(cov(data)))

data.princomp <- princomp(data, cor=T)
plot(data.princomp)

print(summary(data.princomp))

l.components <- 1:4
axes.selected <- data.princomp$scores[,l.components]
plot(axes.selected[,1:2])
plot3d(axes.selected[, c(1,2,3)], main=paste("Component", c(1,2,3), collapse = " "))
plot3d(axes.selected[, c(1,2,4)], main=paste("Component", c(1,2,4), collapse = " "))
plot3d(axes.selected[, c(1,3,4)], main=paste("Component", c(1,3,4), collapse = " "))
plot3d(axes.selected[, c(2,3,4)], main=paste("Component", c(2,3,4), collapse = " "))


# Euclidean distance

d.euclidean <- dist(data, method="euclidean")

clust.hierarchy.single.e <- hclust(d.euclidean, method="single")
plot(clust.hierarchy.single.e, labels=F)
print(paste("Cophenetic cor:", cor(cophenetic(clust.hierarchy.single.e), 
                                   d.euclidean)))

clust.hierarchy.complete.e <- hclust(d.euclidean, method="complete")
plot(clust.hierarchy.complete.e, labels=F)
print(paste("Cophenetic cor:", cor(cophenetic(clust.hierarchy.complete.e), 
                                   d.euclidean)))

mfrow3d(1, 3)
complete.e.2 <- cutree(clust.hierarchy.complete.e, k=2)
plot3d(axes.selected[, c(1,2,3)], main=paste("Component", c(1,2,3), collapse = " "),
       col=c('red', 'green')[complete.e.2])
plot3d(axes.selected[, c(1,2,4)], main=paste("Component", c(1,2,4), collapse = " "),
       col=c('red', 'green')[complete.e.2])

complete.e.3 <- cutree(clust.hierarchy.complete.e, k=3)
plot3d(axes.selected[, c(1,2,3)], main=paste("Component", c(1,2,4), collapse = " "),
       col=c('red', 'green', 'blue')[complete.e.3])

complete.e.5 <- cutree(clust.hierarchy.complete.e, k=5)
plot3d(axes.selected[, c(1,3,4)], main=paste("Component", c(1,3,4), collapse = " "),
       col=c('red', 'green', 'blue', 'orange', 'purple')[complete.e.5])
mfrow3d(1, 1)

clust.hierarchy.average.e <- hclust(d.euclidean, method="average")
plot(clust.hierarchy.average.e, labels=F)
print(paste("Cophenetic cor:", cor(cophenetic(clust.hierarchy.average.e), 
                                   d.euclidean)))

mfrow3d(1, 3)
average.e.2 <- cutree(clust.hierarchy.average.e, k=2)
plot3d(axes.selected[, c(1,2,3)], main=paste("Component", c(1,2,3), collapse = " "),
       col=c('red', 'green')[average.e.2])
plot3d(axes.selected[, c(1,2,4)], main=paste("Component", c(1,2,4), collapse = " "),
       col=c('red', 'green')[average.e.2])

average.e.3 <- cutree(clust.hierarchy.average.e, k=3)
plot3d(axes.selected[, c(1,2,3)], main=paste("Component", c(1,2,4), collapse = " "),
       col=c('red', 'green', 'blue')[average.e.3])

average.e.4 <- cutree(clust.hierarchy.average.e, k=4)
plot3d(axes.selected[, c(1,3,4)], main=paste("Component", c(1,3,4), collapse = " "),
       col=c('red', 'green', 'blue', 'orange', 'purple')[average.e.4])
mfrow3d(1, 1)

# x, y -- vectors from sample
# s.vect -- sample variances by factor
mahalanobis.dist <- function(x, y, s.vect)
{
  delta.sq <- data.frame((x - y)^2)
  sqrt(apply(delta.sq / s.vect, 1, sum))
}

mahalanobis.dist.matr <- function(x.matr)
{
  s.x <- apply(x.matr, 2, var)
  m.d <- function(i, j) { mahalanobis.dist(x.matr[i,], x.matr[j,], s.x) }
  idx <- 1:nrow(x.matr)
  d.result <- as.dist(outer(idx, idx, m.d))
  d.result
}

# Mahalanobis distance

d.mahalanobis <- mahalanobis.dist.matr(data)

clust.hierarchy.single.m <- hclust(d.mahalanobis, method="single")
plot(clust.hierarchy.single.m, labels=F)
print(paste("Cophenetic cor:", cor(cophenetic(clust.hierarchy.single.m), 
                                   d.mahalanobis)))

clust.hierarchy.complete.m <- hclust(d.mahalanobis, method="complete")
plot(clust.hierarchy.complete.m, labels=F)
print(paste("Cophenetic cor:", cor(cophenetic(clust.hierarchy.complete.m), 
                                   d.mahalanobis)))

####
mfrow3d(1, 2)
complete.m.2 <- cutree(clust.hierarchy.complete.m, k=2)
plot3d(axes.selected[, c(1,2,3)], main=paste("Component", c(1,2,3), collapse = " "),
       col=c('red', 'green')[complete.m.2])
plot3d(axes.selected[, c(1,2,4)], main=paste("Component", c(1,2,4), collapse = " "),
       col=c('red', 'green')[complete.m.2])

complete.m.4 <- cutree(clust.hierarchy.complete.m, k=4)
plot3d(axes.selected[, c(1,3,4)], main=paste("Component", c(1,3,4), collapse = " "),
       col=c('red', 'green', 'blue', 'orange', 'purple')[complete.m.4])
mfrow3d(1, 1)
####

clust.hierarchy.average.m <- hclust(d.mahalanobis, method="average")
plot(clust.hierarchy.average.m, labels=F)
print(paste("Cophenetic cor:", cor(cophenetic(clust.hierarchy.average.m), 
                                   d.mahalanobis)))

clust.hierarchy.average.m.4 <- cutree(clust.hierarchy.average.m, k=4)
plot3d(axes.selected[, c(1,3,4)], main=paste("Component", c(1,3,4), collapse = " "),
       col=c('red', 'green', 'blue', 'orange', 'purple')[clust.hierarchy.average.m.4])

# Spectral

mfrow3d(1,2)

clear3d()
clust.spectral.2 <- specc(data, centers=2)
plot3d(axes.selected[, c(1,2,3)], main="On raw",
       col=c("red", "blue")[clust.spectral.2])
clust.spectral.2.s <- specc(axes.selected[,c(1,2,3)], centers=2)
plot3d(axes.selected[, c(1,2,3)], main="On spectral",
       col=c("red", "blue")[clust.spectral.2.s])

clear3d()
clust.spectral.4 <- specc(data, centers=4)
plot3d(axes.selected[, c(1,3,4)], main="On raw",
       col=c("red", "blue", "green", "orange")[clust.spectral.4])
clust.spectral.4.s <- specc(axes.selected[,1:2], centers=4)
plot3d(axes.selected[, c(1,3,4)], main="On spectral",
       col=c("red", "blue", "green", "orange")[clust.spectral.4.s])

mfrow3d(1,1)
