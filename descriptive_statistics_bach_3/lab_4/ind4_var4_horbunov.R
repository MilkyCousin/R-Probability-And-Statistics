library(corrplot)
library(qgraph)

workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

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

df.clo.names <- c(
  'flir','flr','fls','fmc','fosl','frx','fslr','fti','ftr','gas',
  'gci','gd','ge','glid','gls','glw','gme','gnw','goog','gpc'
)

df.clo.full <- form_dataset('./datasets', 6)

df.lgr.full = data.frame(df.clo.full)

corr.clo.pearson <- cor(
  df.clo.full[,-1], method = 'pearson'
)

corrplot.mixed(corr.clo.pearson, title='clo pearson')

corrplot(
  corr.clo.pearson,
  order = 'hclust',
  addrect=4, title='clo pearson, highlighted'
)

pairs(
  corr.clo.pearson[
    ,c('fmc','fti','glw','goog')
  ]
)

pairs(
  corr.clo.pearson[
    ,c('goog', 'flir', 'gd')
    ]
)

corr.clo.spearman <- cor(
  df.clo.full[,-1], method = 'spearman'
)

corrplot(
  corr.clo.spearman,
  order = 'hclust',
  addrect = 5, title='clo spearman, highlighted'
)

corr.clo.kendall <- cor(
  df.clo.full[,-1], method = 'kendall'
)

corrplot.mixed(corr.clo.kendall, title='clo kendall')

corrplot(
  corr.clo.kendall,
  order = 'hclust',
  addrect=5, title='clo kendall, highlighted'
)

pairs(
  corr.clo.kendall[
    ,c('fmc','fti','glw','goog')
    ]
)

pairs(
  corr.clo.kendall[
    ,c('goog', 'flir', 'gd')
    ]
)

qgraph(
  corr.clo.kendall,
  layout='spring', 
  labels=colnames(corr.clo.kendall)
)

# log-return

df.log.ret <- log(df.lgr.full[-nrow(df.lgr.full),]/df.lgr.full[-1,])

corr.log.pearson <- cor(
  df.log.ret, method='pearson'
)

corr.log.kendall <- cor(
  df.log.ret, method='kendall'
)

corr.log.spearman <- cor(
  df.log.ret, method='spearman'
)

corrplot.mixed(corr.log.pearson, title='log pearson')
corrplot.mixed(corr.log.kendall, title='log kendall')
corrplot.mixed(corr.log.spearman, title='log spearman')

corrplot(
  corr.log.spearman,
  order='hclust', addrect=6, title='log spearman, highlighted'
)