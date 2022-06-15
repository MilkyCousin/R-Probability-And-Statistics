set.seed(0)

t.optim <- NULL

C.vect <- NULL

mu.matr <- NULL
d.matr <- NULL
c.matr <- NULL

mu <- function(x.ts, s, t)
{
  mean(x.ts[s:t])
}

d <- function(x.ts, s, t)
{
  sum((x.ts[s:t] - mean(x.ts[s:t]))^2)
}

recursive.mu <- function(x.ts, mu.st, s, t)
{
  ((t - s + 1)*mu.st + x.ts[t+1])/(t - s + 2)
}

recursive.d <- function(x.ts, d.st, mu.st, mu.st1, s, t)
{
  d.st + (t - s + 1) * (mu.st - mu.st1)^2 + (x.ts[t+1] - mu.st1)^2
}

nested.loop <- function(x.ts, K, T.stamp, t.j=c(1))
{
  k <- length(t.j)
  if(k < K)
  {
    if(k == 1)
    {
      for(t.k in (t.j[k]):(T.stamp - K + k))
      {
        C.vect[1] <<- d.matr[t.k,1]
        nested.loop(x.ts, K, T.stamp, c(t.j,t.k+1))
      }
    }
    else
    {
      for(t.k in (t.j[k]):(T.stamp - K + k))
      {
        C.vect[k] <<- C.vect[k-1] + d.matr[t.k,t.j[k]]
        #if((K > 3) & (t.k < 4)) print(c(c.matr[t.k,k], C.vect[k]))
        if(c.matr[t.k,k] >= C.vect[k])
        {
          c.matr[t.k,k] <<- C.vect[k]
          nested.loop(x.ts, K, T.stamp, c(t.j,t.k+1))
        }
      }
    }
  }
  else
  {
    C.vect[k] <<- C.vect[k-1] + d.matr[T.stamp,t.j[k]]
    if(c.matr[T.stamp,k] >= C.vect[k])
    {
      c.matr[T.stamp,k] <<- C.vect[k]
      t.optim <<- c(t.j,T.stamp)
    }
  }
}

intermediate.code <- function(x.ts, K.max = length(x.ts)-1)
{
  T.stamp <- length(x.ts)
  t.opt <- list()
  # row - t, col - k
  
  mu.matr <<- matrix(nrow=T.stamp,ncol=T.stamp)
  for(k in 1:T.stamp)
  {
    mu.matr[k,k] <<- x.ts[k] #mu(x.ts, k, k)
  }
  for(k in 1:(T.stamp-1))
  {
    for(t in (k+1):T.stamp)
    {
      #print(recursive.mu(x.ts, mu.matr[t-1,k], k, t-1))
      mu.matr[t,k] <<- recursive.mu(x.ts, mu.matr[t-1,k], k, t-1)
    }
  }
  
  print("means were calculated.")
  
  d.matr <<- matrix(nrow=T.stamp,ncol=T.stamp)
  for(k in 1:T.stamp)
  {
    d.matr[k,k] <<- 0 # (x.ts[k] - mu(x.ts, k, k))^2 = 0
  }
  for(k in 1:(T.stamp-1))
  {
    for(t in (k+1):T.stamp)
    {
      d.matr[t,k] <<- recursive.d(x.ts, d.matr[t-1,k], mu.matr[t-1,k], mu.matr[t,k], k, t-1)
    }
  }
  
  print("errors were calculated.")
  c.matr <<- d.matr
  
  for(K in 2:K.max)
  {
    C.vect <<- matrix(nrow=1,ncol=K)
    nested.loop(x.ts, K, T.stamp)
    t.opt[[K-1]] <- t.optim
    
    t.optim <<- NULL
  }
  t.opt
}

test.0 <- function()
{
  T.full <- 10000
  x.test <- c(rep(1, T.full * 0.25), rep(1.5, T.full * 0.5), rep(0.5, T.full * 0.25)) + rnorm(T.full, sd=0.25)
  
  t0 <- Sys.time()
  t.matr <- intermediate.code(x.test, K.max = 4)
  t1 <- Sys.time()
  
  print(t1-t0)
  
  for(u in 1:length(t.matr))
  {
    plot(x.test, type='l'); grid()
    results <- t.matr[[u]]
    for(i in 1:(length(results)-1)) 
    {
      #print(results)
      segments(results[i],mu(x.test,results[i],results[i+1]), results[i+1], mu(x.test,results[i],results[i+1]), lwd=2, col='darkred')
      #print(mu(x.test,results[i],results[i+1]))
    }
  }
}

test.1 <- function()
{
  m <- c(4, 5, 4.4, 5, 4.3, 3.6, 4.4, 3.1, 3.5, 4.9)
  x.test <- c()
  N <- 1000
  for(m.i in m)
  {
    x.test <- c(x.test, rep(m.i, N) + rnorm(N))
  }
  
  t0 <- Sys.time()
  t.matr <- intermediate.code(x.test, K.max = 10)
  t1 <- Sys.time()
  
  print(t1-t0)
  
  for(u in 1:length(t.matr))
  {
    plot(x.test, type='l'); grid()
    results <- t.matr[[u]]
    for(i in 1:(length(results)-1)) 
    {
      #print(results)
      segments(results[i],mu(x.test,results[i],results[i+1]), results[i+1], mu(x.test,results[i],results[i+1]), lwd=2, col='darkred')
      #print(mu(x.test,results[i],results[i+1]))
    }
  }
}