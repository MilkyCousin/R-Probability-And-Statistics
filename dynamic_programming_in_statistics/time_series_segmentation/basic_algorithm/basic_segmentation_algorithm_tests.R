# В этом файле описано тестирование работы алгоритма сегментации на некоторых
# данных, созданных вручную. Три выборки, как было указано в тексте моей работы.
# Не совсем красиво оформлено, скорость работы программы приемлемая.  

workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()
# Для полноценной работы необходимо хранить оба файла в одной директории
if(!exists("segment.err", mode = "function")) source("basic_segmentation_algorithm_horbunov_iii_cstat.R")

# utility
# later: better plotting function

plot.two.samples.as.tss <- function(x, y, gp = list(), lg = c())
{
  ts.plot(ts(x, start = 1, end = length(x)),
          ts(y, start = 1, end = length(y)),
          gpars = gp)
  legend("topleft", legend = lg, col = c("red", "black"), lty = 1)
}

step.wise.data <- function(start, end, m, step = 1)
{
  order <- seq(start, end, step)
  r <- c()
  for(k in order)
  {
    r <- c(r, rep(k, m))
  }
  r
}

plot.as.ts <- function(x.sample, gpars.list = list())
{
  ts.x.s <- ts(x.sample, start = 1, end = length(x.sample))
  ts.plot(ts.x.s, gpars = gpars.list)
}

# main testing function, based on methods written in my work

tests.basic.full <- function(x, K.max, N.val, N, n.v.exp, v.v.sdev, boolean.n = F, boolean.i = F, boolean.d = F)
{
 l.sdev <- length(v.v.sdev)
 
 r.v.p.z <- vector(mode = "numeric", length = l.sdev)
 r.v.b.z <- vector(mode = "numeric", length = l.sdev)
 
 ###################
 
 for(i in 1:l.sdev)
 {
   print(sprintf("%d. current sdev: %g", i, v.v.sdev[i]))
   
   # plots sample and its modified version
   plot.separated(x, n.v.exp, v.v.sdev[i]) # ? enormous values
   
   r.c <- test.basic.bmetric.and.proportions(x, K.max, N.val, N, n.v.exp, v.v.sdev[i], boolean.n, boolean.i)
   
   r.v.p.z[i] <- r.c[[1]]
   r.v.b.z[i] <- r.c[[2]]
 }
 
 ###################
  
 plot(v.v.sdev, r.v.p.z, type = "l", main = "proportion loss -- sdev")
 plot(v.v.sdev, r.v.b.z, type = "l", main = "beeferman metric value -- sdev")
 
 matplot(v.v.sdev, cbind(r.v.p.z, r.v.b.z), col = c("red", "blue"), type = "l", ylab = "errors", xlab = "sdev")
 legend("topleft", legend = c("proportion loss", "beeferman metric"), col = c("red", "blue"), lty = 1)
}

tests.basic.full.sdev.num <- function(x, K.max, N.val, N, n.v.exp, n.v.sdev, boolean.n = F, boolean.i = F)
{
  plot.separated(x, n.v.exp, n.v.sdev)    #(x, K.max, N.val, N, n.v.exp, v.v.sdev[i], boolean.n, boolean.i)
  unlist(test.basic.bmetric.and.proportions(x, K.max, N.val, N, n.v.exp, n.v.sdev, boolean.n, boolean.i))
}


t.0 <- function(x, k.max, N.val, N, m.n, sdev.v, b.n = F, b.m.i = F)
{
  tests.basic.full.sdev.num(x, k.max, N.val, N, m.n, sdev.v, b.n, b.m.i)
}

t.1 <- function(N.val = 1000, N = 1, mul = 1, m.n = 0, sdev.v = 1, b.n = F, b.m.i = F)
{
  z <- c(rep(0.4, 35 * mul), rep(0.6, 30 * mul), rep(0.8, 35 * mul))
  z.kmax <- 3
  t.0(z, z.kmax, N.val, N, m.n, sdev.v, b.n, b.m.i)
}

t.2 <- function(N.val = 1000, N = 1, mul = 1, m.n = 0, sdev.v = 1, b.n = F, b.m.i = F)
{
  z <- c(rep(0.25, 25 * mul), rep(0.65, 50 * mul), rep(0.45, 25 * mul))
  z.kmax <- 3
  t.0(z, z.kmax, N.val, N, m.n, sdev.v, b.n, b.m.i)
}

t.3 <- function(N.val = 1000, N = 1, mul = 1, m.n = 0, sdev.v = 1, b.n = F, b.m.i = F)
{
  z <- c(step.wise.data(0.1, 0.4, 15 * mul, 0.15), rep(0.6, 10 * mul), step.wise.data(0.4, 0.1, 15 * mul, -0.15))
  z.kmax <- 7
  t.0(z, z.kmax, N.val, N, m.n, sdev.v, b.n, b.m.i)
}

### testing
# note: up to 10K elements
# suppose that variable 'mul' is a natural number

plot.separated <- function(x, m.n, sdev.v)
{
  x.n <- modify.by.rnorm(x, m.n, sdev.v)
  plot.two.samples.as.tss(x, x.n, gp = list(col = c("red", "black")), lg = c("μ = 0, σ = 0", sprintf("μ = %g, σ = %g", m.n, sdev.v)))
}

# arguments:
# x - sample | time series
# k.max - max number of parts
# N - number of experiments
# m.n - expected value for rnorm()
# sdev.v - vector of sdev's for rnorm()
# b.m.i - whether to use time efficient beeferman metric function or not

main_t <- function(x, k.max, N.val, N, m.n, sdev.v, b.m.i = F)
{
  # main, normed
  tests.basic.full(x, k.max, N.val, N, m.n, sdev.v, T, b.m.i)
}

main_1 <- function(N = 1000, block.len = 1, mul = 1, b.m.i = F)
{
  # A - sample test
  z <- c(rep(0.4, 35 * mul), rep(0.6, 30 * mul), rep(0.8, 35 * mul))
  z.kmax <- 3
  z.m = 0
  sdev.z <- seq(0, 0.25, 0.05)
  main_t(z, z.kmax, block.len, N, z.m, sdev.z, b.m.i)
}

main_2 <- function(N = 1000, block.len = 1, mul = 1, b.m.i = F)
{
  # B - sample test
  z <- c(rep(0.25, 25 * mul), rep(0.65, 50 * mul), rep(0.45, 25 * mul))
  z.kmax <- 3
  z.m = 0
  sdev.z <- seq(0, 0.25, 0.05)
  main_t(z, z.kmax, block.len, N, z.m, sdev.z, b.m.i)
}

main_3 <- function(N = 1000, block.len = 1, mul = 1, b.m.i = F)
{
  # C - sample test
  z <- c(step.wise.data(0.1, 0.4, 15 * mul, 0.15), rep(0.6, 10 * mul), step.wise.data(0.4, 0.1, 15 * mul, -0.15))
  z.kmax <- 7
  z.m = 0
  sdev.z <- seq(0, 0.25, 0.05)
  main_t(z, z.kmax, block.len, N, z.m, sdev.z, b.m.i)
}
