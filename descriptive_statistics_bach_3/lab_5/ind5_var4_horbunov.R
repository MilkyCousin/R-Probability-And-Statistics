library(rgl)

workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

lab5.operations <- function(
  numeric.data, 
  boolean.enable.prints=F, boolean.enable.cor=F, boolean.3d = T
)
{
  pairs(
    numeric.data, cex = 0.1
  )
  
  m.comp.cov <- princomp(
    numeric.data,
    cor = boolean.enable.cor
  )
  
  if(boolean.enable.prints)
  {
    print('sdev:')
    print(
      m.comp.cov$sdev
    )
    print('loadings:')
    print(
      m.comp.cov$loadings
    )
    print('scores:')
    print(
      m.comp.cov$scores
    )
  }
  
  m.comp.summary <- summary(m.comp.cov)
  print(m.comp.summary)
  
  if(boolean.3d)
  {
    boolean.report <- T
    if(boolean.report)
    {
      i = 3; j = 7
    }
    else
    {
      i <- as.numeric(readline(prompt = "ind1="))
      j <- as.numeric(readline(prompt = "ind2="))
    }
    
    plot(m.comp.cov$scores[,i:j])
    plot3d(m.comp.cov$scores[,c(3,7)])
  }
}

text.data <- read.table('./multi/F4p.txt', header = F)

lab5.operations(text.data)

