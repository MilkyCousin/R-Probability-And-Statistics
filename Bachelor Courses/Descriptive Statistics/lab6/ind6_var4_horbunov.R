library(rgl)

form_dataset <- function(path.to.files, ind)
{
  list.of.paths <- list.files(
    path.to.files,
    full.names = T
  )
  
  print(
    list.of.paths
  )
  
  extracted.data <- lapply(
    list.of.paths, 
    function(t)
    {
      t.raw <- read.csv(
        file = t, header = F
      )[,c(1,ind)];
      colnames(t.raw) <- c(
        'data', unlist(strsplit(t, '[_.]'))[3]
      );
      t.raw
    }
  )
  
  result.frame <- Reduce(
    function(i, j)
    {
      merge(i, j, by='data')
    },
    extracted.data
  )
  
  result.frame
}

workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

set.distance <- function(n.data, s.method)
{
  r.dist <- dist(n.data, method=s.method)
  r.dist
}

set.hclust <- function(n.dist, s.method)
{
  r.hclust <- hclust(n.dist, method=s.method)
  r.hclust
}

set.plot <- function(n.data, s.title, boolean.as.dendrogram = F)
{
  plot(
    if(boolean.as.dendrogram) as.dendrogram(n.data) else n.data, 
    main=s.title, leaflab="none"
  )
}
# 
# lab6.operations <- function(numeric.data, boolean.dendrogram = F)
# {
#   m.comp <- princomp(numeric.data)
#   m.comp.summary <- summary(m.comp)
#   
#   plot(numeric.data, cex=0.1)
#   print(m.comp.summary)
#   
#   get.dendrogram <- function(m.comp, i, j, q, k.value, s.title='')
#   {
#     m.scores.slice <- m.comp$scores[,i:j]
#     m.comp.dist.euc <- set.distance(m.scores.slice, 'euclidean')
#     m.comp.dist.max <- set.distance(m.scores.slice, 'maximum')
#     m.comp.dist.man <- set.distance(m.scores.slice, 'manhattan')
#     
#     s.methods <- c(
#       'single', 'complete', 'average'
#     )
#     
#     m.hclust.euc <- set.hclust(m.comp.dist.euc, s.methods[q])
#     m.hclust.max <- set.hclust(m.comp.dist.max, s.methods[q])
#     m.hclust.man <- set.hclust(m.comp.dist.man, s.methods[q])
#     
#     set.plot(m.hclust.euc, "euclidean dist hclust", boolean.dendrogram)
#     set.plot(m.hclust.euc, "maximum dist hclust", boolean.dendrogram)
#     set.plot(m.hclust.euc, "l_1 dist hclust", boolean.dendrogram)
#     
#     m.hclust.groups <- cutree(m.hclust.euc, k=k.value)
#     
#     s.colors <- c(
#       'green', 'red', 'blue', 'purple'
#     )
#     
#     print(s.colors[1:k.value])
#     
#     plot(
#       m.comp$scores[,c(i,j)], cex=0.5,
#       col=s.colors[1:k.value][m.hclust.groups],
#       main=s.title
#     )
#     
#     plot3d(m.comp$scores[,c(i,j)], 
#            col=s.colors[1:k.value][m.hclust.groups],size=4)
#   }
#   
#   list.indices <- list(
#     c(5, 7, 2, 4)
#   )
#   
#   for(pair in list.indices)
#   {
#     get.dendrogram(
#       m.comp, pair[1], pair[2], pair[3], pair[4],
#       paste(pair[1], pair[2], sep=":")
#     )
#   }
# }
#   
# text.data <- read.table('./multi/F4p.txt', header = F)
# 
# lab6.operations(text.data, boolean.dendrogram = T)

 colors.global <- c('red', 'green', 'blue', 'purple', 'orange', 'pink', 'gray', 'yellow')
#
 df.clo.values <- form_dataset('./datasets', 6)
 df.clo.m.comp <- princomp(df.clo.values[,-1])
#
 i.1 = 3; j.1 = 4

# corrmat <- cor(df.clo.values, method='kendall')
# corrplot(corrmat, order='hclust')
#
# plot(
#   df.clo.m.comp$scores[,c(i.1,j.1)]
# )


 df.clo.dist <- dist(df.clo.m.comp$scores[,i.1:j.1])
 df.clo.hclust <- hclust(df.clo.dist, method="single")
 plot(as.dendrogram(df.clo.hclust), leaflab='none', main='clo hclust')

 df.clo.hclust.groups <- cutree(df.clo.hclust, k=4)
 plot(
   df.clo.m.comp$scores[,c(i.1,j.1)],
   col=colors.global[df.clo.hclust.groups],
   main='clust df clo'
 )
 plot3d(
   df.clo.m.comp$scores[,c(i.1,j.1)],
   col=colors.global[df.clo.hclust.groups],size=4
)