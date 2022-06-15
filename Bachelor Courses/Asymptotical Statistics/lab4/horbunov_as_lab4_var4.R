workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

dat <- read.table("D4.txt", header = T)
x <- as.double(dat[,1])
y <- as.double(dat[,2])

plot(x, y)

T <- 3

nls.xy <- nls(y ~ cos(b * x) * log(a + x), 
              start = list(
                a = 1,
                b = 2*pi/T
              ),
              lower = list(
                a = 10e-5,
                b = 10e-5
              ),
              algorithm = "port")
print(summary(nls.xy))

pars <- nls.xy$m$getPars()
a <- pars[1]
b <- pars[2]

plot(x, y)

t <- seq(0.01,10,0.01)
lines(t, cos(b * t) * log(t + a), col = 'red', lwd = 2)

u <- nls.xy$m$resid()