library(gridExtra)
library(latticeExtra)
library(viridisLite)

na.erase <- function(x)
{
  Filter(function(t){return(!is.na(t))}, x)
}

dist.matr <- function(x.sample, y.sample, dist.method = 'euclidean')
{
  N <- length(x.sample)
  M <- length(y.sample)

  cost.local <- matrix(nrow = N, ncol = M)
  
  for(i in 1:N)
  {
    for(j in 1:M)
    {
      cost.local[i,j] = dist(c(x.sample[i],y.sample[j]), method = dist.method)
    }
  }
  
  cost.local
}

accumulated.cost.matrix <- function(x.sample, y.sample, c.local)
{
  N <- length(x.sample)
  M <- length(y.sample)
  
  dtw.matrix <- matrix(nrow = N, ncol = M)
  
  dtw.matrix[1,1] <- 0
  for(i in 2:N)
  {
    dtw.matrix[i,1] <- dtw.matrix[i-1,1] + c.local[i,1]
  }
  for(j in 2:M)
  {
    dtw.matrix[1,j] <- dtw.matrix[1,j-1] + c.local[1,j]
  }
  for(i in 2:N)
  {
    for(j in 2:M)
    {
      dtw.matrix[i,j] = c.local[i,j] + min(dtw.matrix[i-1,j], dtw.matrix[i,j-1], dtw.matrix[i-1,j-1])
    }
  }
  
  dtw.matrix
}

optimal.warping.path <- function(dtw.m)
{
  path.dtw.x <- numeric(0)
  path.dtw.y <- numeric(0)
  
  i <- nrow(dtw.m)
  j <- ncol(dtw.m)
  
  path.dtw.x <- c(path.dtw.x, i)
  path.dtw.y <- c(path.dtw.y, j)
  
  while((i>1) && (j>1))
  {
    if(i==1)
    {
      j <- j-1
    }
    else if(j==1)
    {
      i <- i-1
    }
    else
    {
      dtw.min <- min(c(dtw.m[i-1,j],dtw.m[i,j-1],dtw.m[i-1,j-1]))
      if(dtw.m[i-1,j] == dtw.min)
      {
        i <- i-1
      }
      else if(dtw.m[i,j-1] == dtw.min)
      {
        j <- j-1
      }
      else
      {
        i <- i-1
        j <- j-1
      }
      path.dtw.x <- c(path.dtw.x, i)
      path.dtw.y <- c(path.dtw.y, j)
    }
  }
  
  path.dtw.x <- c(path.dtw.x, 1)
  path.dtw.y <- c(path.dtw.y, 1)

  data.frame(path.x = path.dtw.x, path.y = path.dtw.y)
}

plt.dtw.path <- function(obj)
{
  print(
    levelplot(
      obj$dist.matr,
    col.regions = viridis(100), 
    cuts = 64,
    main = paste('DTW warping path plot, metric:', obj$dist.method),
    xlab = 'Query index',
    ylab = 'Template index') + contourplot(obj$dist.matr) + xyplot(
      obj$opt.path[,2] ~ obj$opt.path[,1], 
      type="l", lwd=1, col='red'
      ) 
  )
}

plt.templated <- function(x, y, t, z, title, legend.v, colors.v, pos.v = "bottomleft")
{
  mn = min(c(x,y,z))
  mx = max(c(x,y,z))
  plot(x, type='l', col=colors.v[1], xlim=c(0,max(length(x),length(y))+1),
       main=title, xlab='index', ylab='value', lwd=1.5, ylim=c(mn,mx)) + lines(
    y, type='l', col=colors.v[2], lwd=1.5) + lines(t, z, col=colors.v[3], lty=2, lwd=1)
  legend(pos.v,legend=legend.v, col=colors.v, lwd=1)
}

plt.dtw.align <- function(obj)
{
  lg.v <- c('x', 'y', 'align')
  cl.v <- c('red', 'blue', 'darkgreen')
    
  plt.templated(obj$query, obj$template, obj$opt.path[,2], obj$query[obj$opt.path[,1]],
                paste('x and y plot after alignment (by x), metric:', obj$dist.method),
                lg.v, cl.v)
  
  plt.templated(obj$query, obj$template, obj$opt.path[,1], obj$template[obj$opt.path[,2]],
                paste('x and y plot after alignment (by y), metric:', obj$dist.method),
                lg.v, cl.v)
}

plt.dtw.mapping <- function(obj)
{
  mn = min(c(obj$query, obj$template))
  mx = max(c(obj$query, obj$template))
  
  plot(
    obj$query, type='l', col='red', lwd=1.5, main='DTW mapping', 
    xlab='index', ylab='value', ylim=c(mn,mx), 
    xlim=c(0,max(length(obj$query),length(obj$template))+1),
    ) + lines(obj$template, col='blue', lwd=1.5)
  
  segments(
    obj$opt.path[,2], obj$template[obj$opt.path[,2]],
    obj$opt.path[,1], obj$query[obj$opt.path[,1]],  # q <- t
    col='dimgray', lty=2)
}

func.calculate <- function(distances.matr, path.aligned)
{
  f.v <- 0
  for(j in 1:nrow(path.aligned))
  {
    f.v <- f.v + distances.matr[path.aligned[j,][[1]],path.aligned[j,][[2]]]
  }
  f.v/nrow(path.aligned)
}

basic.dtw <- function(x.series, y.series, method.dist='euclidean')
{
  distances <- dist.matr(x.series, y.series, method.dist)
  accumulated <- accumulated.cost.matrix(x.series, y.series, distances)
  path.optimal <- optimal.warping.path(accumulated)
  
  dtw.results <- list(
    query = x.series,
    template = y.series,
    dist.method = method.dist,
    dist.matr = distances, 
    acc.matr = accumulated, 
    opt.path = path.optimal,
    opt.func = func.calculate(distances, path.optimal)
    )
  
  class(dtw.results) <- "dtwmodel"
  dtw.results
}

plot.dtwmodel <- function(obj, ...)
{
  plt.dtw.align(obj)
  plt.dtw.mapping(obj)
  plt.dtw.path(obj)
}

#i <- seq(0, 2*3.14, len=100);
#x <- cos(i) + runif(100) / 10;
#y <- sin(i)

#r <- basic.dtw(x,y)