library(raster)

workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

func <- function()
{
  brazil <- getData(
    "GADM", country="BRA", level=1
  )
  
  save(
    brazil, file="bramap.Rdata"
  )
}

load("bramap.Rdata")
print(
  brazil$NAME_1
)
write.csv(
  brazil$NAME_1, file="brazil_states.csv"
)
