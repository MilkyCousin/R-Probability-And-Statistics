library('kernlab')
library('rgl')

options(rgl.printRglwidget = TRUE)

data <- read.table("F4p.txt", header=F)
plot(data, cex=0.25)

print(diag(cov(data)))

data.princomp <- princomp(data, cor=T)
plot(data.princomp)

print(summary(data.princomp))

axes.selected <- data.princomp$scores
pairs(axes.selected, cex=0.25)

selected.1 <- axes.selected[,c(1,3,4)]
plot3d(selected.1)

mfrow3d(2,2)
plot3d(axes.selected[,c(3,6,10)])
plot3d(axes.selected[,c(3,6,10)])
plot3d(axes.selected[,c(3,6,10)])
plot3d(axes.selected[,c(3,6,10)])
mfrow3d(1,1)

kmeans.spec <- kmeans(axes.selected[,c(3,6,10)], centers = 4)
print(kmeans.spec)

library('factoextra')
fviz_nbclust(
  x=axes.selected[,c(3,6,10)],
  FUNcluster=kmeans,
  k.max=10,
  method="silhouette"
)
fviz_nbclust(
  x=axes.selected[,c(3,6,10)],
  FUNcluster=kmeans,
  k.max=10,
  method="wss"
)
pal <- c('red', 'blue', 'green', 'orange')

km.2 <- kmeans(x=axes.selected[,c(3,6,10)],centers=2)
pairs(axes.selected, cex=0.25, col=pal[km.2$cluster])

km.3 <- kmeans(x=axes.selected[,c(3,6,10)],centers=3)
pairs(axes.selected, cex=0.25, col=pal[km.3$cluster])

km.4 <- kmeans(x=axes.selected[,c(3,6,10)],centers=4)
pairs(axes.selected, cex=0.25, col=pal[km.4$cluster])

library('fossil')

set.seed(0)
specc.2 <- specc(x=axes.selected[,c(3,6,10)], centers=2)
specc.3 <- specc(x=axes.selected[,c(3,6,10)], centers=3)
specc.4 <- specc(x=axes.selected[,c(3,6,10)], centers=4)



pairs(axes.selected, col=pal[specc.4], cex=0.25)

plot3d(selected.1, col=pal[specc.4])
