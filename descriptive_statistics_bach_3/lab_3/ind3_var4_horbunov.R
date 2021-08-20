library("sp")
library("maptools")

#todo: chk

workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

bra_data <- read.csv(
  "geoMap.csv", header=T,
  stringsAsFactors=F
)

  bra_names <- read.csv(
  "brazil_states.csv", header=T
)

full_data <- merge(
  bra_names, bra_data,
  by="Region", all=T
)

map_print <- function(ratio_1, string_v, str_vect)
{
  numcol <- 7
  brkMP <- c(
    min(ratio_1)-0.001,
    seq(
      min(ratio_1)+0.001,
      max(ratio_1)+0.001,
      length.out=numcol+1
    )
  )
  intMP <- numcol + 2 - findInterval(ratio_1, brkMP)
  
  palette(
    c(heat.colors(numcol, alpha = 1), rgb(1,1,1))
  )
  par(xpd=T)
  par(mar=c(0,0,0,0))
  
  load("bramap.Rdata")
  
  print(string_v)
  print(brkMP)
  
  plot(brazil, col=intMP[order(full_data$x)])
  
  legend(
    "bottomright", title=string_v,
    legend = str_vect,
    fill=c(
      heat.colors(
        numcol, alpha=1
      ),
      rgb(1,1,1)
    ), inset=0, horiz=F
  )
}

ratio_n1 = full_data$macaw/(
  full_data$macaw + full_data$parakeet + full_data$pigeons
)

ratio_n2 = full_data$parakeet/(
  full_data$macaw + full_data$parakeet + full_data$pigeons
)

ratio_n3 = full_data$pigeons/(
  full_data$macaw + full_data$parakeet + full_data$pigeons
)

map_print(ratio_n1, 'macaw %', c('<20%','<25%','<35%','<45%','>53%'))

map_print(ratio_n2, 'parakeet %', c('<11%','<20%','<35%','<45%','>53%'))

map_print(ratio_n3, 'pigeons %', c('<23%','<30%','<40%','<45%','>52%'))
