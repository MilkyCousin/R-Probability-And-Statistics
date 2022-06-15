# В этом файле описан базовый алгоритм сегментации выборки по частям.
# К тому же, присутствует программная реализация оценок (за отклонениями итоговых пропорций; за метрикой Биффермана) работы
# необходимого алгоритма.

#### UTILITY ####

na.erase <- function(x)
{
  Filter(function(t){return(!is.na(t))}, x)
}

### i_e <- c(i1,i2,...in)
### r <- c(i1-is0, i2-i1, ..., in-in-1)

#### DP segmentation algorithm based on segment mean values ####

# attempt to implement dynamic programming segmentation algorithm

# segmentation error computation algorithm, based on segment-mean
segment.err <- function(x.sample, s, t)
{
  x.slice <- x.sample[s:t]
  x.slice.mean <- mean(x.slice) # med?
  result <- sum((x.slice - x.slice.mean) ^ 2)
  result
}

seq.mod <- function(from, to, by)
{
  seq.v <- seq(from, to, by)
  if(tail(seq.v, 1) != to)
  {
    return(c(seq.v, to))
  }
  return(seq.v)
}

segment.optimization.process <- function(x.sample, K.max, N = 1)
{
  T.ind <- length(x.sample)
  s.ind <- seq.mod(1, T.ind, N)
  l.s.ind <- length(s.ind)
  d.err <- matrix(nrow = l.s.ind, ncol = l.s.ind)
  for(t in 1:l.s.ind)
  {
    for(s in 1:t)
    {
      d.err[s, t] <- segment.err(x.sample, s.ind[s], s.ind[t])
    }
  }
  c.matr <- matrix(nrow=l.s.ind, ncol=K.max)
  z.matr <- matrix(nrow=l.s.ind, ncol=K.max)
  for(t in 1:l.s.ind)
  {
    c.matr[t,1] <- d.err[1, t]  # ci(1) = d(1,t)
    z.matr[t,1] = 1
  }
  # optimization part
  for(k in 2:K.max)
  {
    c.matr[1,k] = 0
    for(s in 2:l.s.ind)
    {
      e.vect <- vector(mode="numeric", length=s-1)
      for(t in 1:s-1)
      {
        e.vect[t] <- c.matr[t,k-1] +  d.err[t+1, s]
      }
      c.matr[s,k] <- min(e.vect)
      z.matr[s,k] <- s.ind[which.min(e.vect)]
    }
  }
  # backtracking part
  t.optimal.ind <- matrix(nrow=K.max,ncol=K.max)
  for(i in 1:K.max)
  {
    t.optimal.ind[i, i] = T.ind
    for(j in rev(1:i))
    {
      t.optimal.ind[j-1, i] = z.matr[match(t.optimal.ind[j, i], s.ind), j] # t(k-1)K
    }
  }
  # obtaining result as matrix with indices 4 each column represents optimal indices for a specific number of segments
  t.optimal.ind
}

#### UTILITY, CURRENTLY USING ####

get.proportions.kmax <- function(t.optimal.ind, K.max)
{
  t.K.ind <- na.erase(t.optimal.ind[,K.max])
  proportions.K <- c(t.K.ind[1], diff(t.K.ind))/max(t.K.ind)
  proportions.K
}

modify.by.rnorm <- function(x.sample, exp = 0, sdev = 1)
{
  x.sample.m <- x.sample + rnorm(length(x.sample), exp, sdev)
  x.sample.m
}

generate.samples.rnorm <- function(x.sample, m, exp = 0, sdev = 1)
{
  list.x.sample.m.copies <- list(length = m)
  for (i in 1:m)
  {
    list.x.sample.m.copies[[i]] <- modify.by.rnorm(x.sample, exp, sdev)
  }
  list.x.sample.m.copies
}

binary.delta <- function(t.x, i, j)
{
  k = 1
  for (l in t.x)
  {
    if(((k <= i) && (i <= l)) && ((k <= j) && (j <= l)))
    {
      return(1)
    }
    k = l
  }
  return(0)
}

#### METRICS ####

beeferman.segmentation.metric <- function(t.y, t.z, boolean.debug = F, n = F)
{
  T <- t.y[length(t.y)]
  div <- if(n) (T)^2 else 1
  output.result <- 0
  for(i in 1:T)
  {
    for(j in 1:T)
    {
      output.result = abs(binary.delta(t.y, i, j) - binary.delta(t.z, i, j)) + output.result
    }
  }
  output.result = output.result / div
  output.result
}

beeferman.segmentation.metric.improved <- function(t.y, t.z, boolean.debug = F, n = F)
{
  if(length(t.y) > 1)
    k <- floor(0.5 * mean(c(diff(t.y), diff(t.z))))
  else
    k <- 0
  T <- t.y[length(t.y)]
  div <- if(n) (T-k-1) else 1
  output.result <- 0
  for(i in 1:T-k-1)
  {
    output.result = abs(binary.delta(t.y, i, i + k + 1) - binary.delta(t.z, i, i + k + 1)) + output.result
  }
  output.result = output.result / div
  output.result
}

#### CURRENTLY USING ####

test.basic.bmetric.and.proportions <- function(z.sample, k.max, n, m, exp = 0, sdev = 1, boolean.n = F, boolean.i = F)
{
  t.z <- segment.optimization.process(z.sample, k.max, n)
  
  z.s.m.l <- generate.samples.rnorm(z.sample, m, exp, sdev)
  t.z.s.m.l <- lapply(z.s.m.l, segment.optimization.process, K.max=k.max, N=n)
    
  #### PROPORTIONS ####
  
  r.d.sum <- 0
  
  for(k in 1:k.max)
  {
    p.t.z.s.m.l <- Map(get.proportions.kmax, t.z.s.m.l, k)
    p.t.z <- get.proportions.kmax(t.z, k)
    
    p.t.z.s.m.l <- lapply(p.t.z.s.m.l, function(a, b){return((a-b)^2)}, b = p.t.z)
    v.proportions <- sqrt(unlist(lapply(p.t.z.s.m.l, sum))/k)
    r.d.sum <- r.d.sum + mean(v.proportions)
  }
  
  #### RAND INDEX ####

  r.b.d.sum <- 0
  
  for(k in 1:k.max)
  {
    t.z.s.m.l.c <- lapply(t.z.s.m.l, function(x){return(na.erase(x[,k]))})
    
    if(boolean.i)
      b.dist.l <- unlist(lapply(t.z.s.m.l.c, beeferman.segmentation.metric.improved, t.y = na.erase(t.z[,k]), n = boolean.n))
    else
      b.dist.l <- unlist(lapply(t.z.s.m.l.c, beeferman.segmentation.metric, t.y = na.erase(t.z[,k]), n = boolean.n))
      
    b.dist.m <- mean(b.dist.l)
    
    r.b.d.sum <- r.b.d.sum + b.dist.m
  }
  
  list(r.d.sum/k.max, r.b.d.sum/k.max)
}

# information

j.functional <- function(x.sample, t.vect)
{
  j.value <- 0
  v.ind <- c(1, na.erase(t.vect))
  for(s in 1:(length(v.ind)-1))
  {
    j.value <- j.value + segment.err(x.sample, v.ind[s], v.ind[s+1])
  }
  j.value
}

bic.segmentation <- function(x.sample, t.optimal.matr, K.max)
{
  K <- 1:K.max
  T.v <- t.optimal.matr[1,1]
  j.vect <- c()
  for(k in K)
  {
    j.vect <- c(j.vect, j.functional(x.sample, t.optimal.matr[,k]))
  }
  j.min <- min(j.vect)
  J.estimate <- T.v*log(j.vect/(T.v-1)) + 2*K*log(T.v)
  print('functional outcomes:')
  print(j.vect)
  print('bic outcomes:')
  print(J.estimate)
  c(min(J.estimate), which.min(J.estimate))
}

fast.test.n <- function(mul=1, blocks=1)
{
  z <- c(rep(0.4, 35 * mul), rep(0.6, 30 * mul), rep(0.8, 35 * mul))
  r <- test.basic.bmetric.and.proportions(z.sample=z, k.max=3, n=blocks, m=1000, exp = 0, sdev = .15, boolean.n = T, boolean.i = T)
  print(r)
}

#t1
#m 1    B 1  0.08665775 0.12232713
#m 2.5  B x  0.06342043 0.09451996