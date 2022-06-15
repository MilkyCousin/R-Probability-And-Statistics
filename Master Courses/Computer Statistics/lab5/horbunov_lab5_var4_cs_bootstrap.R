# Фіксуємо зернину
set.seed(1)

g <- function(t)
{
  sqrt(0.5)*abs(t-1)
}

a <- -2
b <- 0.5

m <- -1
s <- sqrt(2)

bootstrap.gls <- function(N = 1000, B = 1000)
{
  b.ols <- c()
  b.gls <- c()
  
  for(i in 1:B)
  {
    X <- data.matrix(rnorm(N, mean=m, sd=s))
    E <- data.matrix(rnorm(N, mean=0, sd=g(X)))
    Y <- a + b*X + E
    
    mX = mean(X)
    mY = mean(Y)
    
    b_ols = cov(X, Y)/var(X)
    a_ols = mY - b_ols * mX
    bo <- data.matrix(c(a_ols, b_ols))
    b.ols <- cbind(b.ols, bo)
    
    X.i <- cbind(1, X)
    
    Yhat_ols <- X.i%*%bo
    U_ols <- Y - Yhat_ols
    
    X.poly <- cbind(1,X,X^2)
    A.poly <- t(X.poly)%*%X.poly
    greeks <- solve(A.poly)%*%t(X.poly)%*%(U_ols^2)
    
    Usq_fitted <- X.poly%*%greeks
    
    W.res <- diag(as.numeric(1 / Usq_fitted))
    W.res.0 <- W.res*((0 < W.res) & (W.res <= 1)) + 1*(W.res > 1)
    
    Aw <- t(X.i)%*%W.res.0%*%X.i
    bw <- solve(Aw)%*%t(X.i)%*%W.res.0%*%Y
    b.gls <- cbind(b.gls, bw)
  }
  
  rownames(b.ols) <- c("a", "b")
  rownames(b.gls) <- c("a", "b")
  
  list(ols.vect = b.ols, gls.vect = b.gls, len = N)
}

operate.lists <- function(bs.list)
{
  ma.ols <- mean(bs.list$ols.vect[1,])
  mb.ols <- mean(bs.list$ols.vect[2,])
  
  ma.gls <- mean(bs.list$gls.vect[1,])
  mb.gls <- mean(bs.list$gls.vect[2,])
  
  bias.a.ols <- a - ma.ols
  bias.b.ols <- b - mb.ols
  
  bias.a.gls <- a - ma.gls
  bias.b.gls <- b - mb.gls
  
  print(paste("Bias, N = ", bs.list$len))
  print("OLS:")
  print(c(bias.a.ols, bias.b.ols))
  print("GLS:")
  print(c(bias.a.gls, bias.b.gls))
  
  sd.a.ols <- sd(bs.list$ols.vect[1,])
  sd.b.ols <- sd(bs.list$ols.vect[2,])
  
  sd.a.gls <- sd(bs.list$gls.vect[1,])
  sd.b.gls <- sd(bs.list$gls.vect[2,])
  
  print(paste("sqrt(Variance), N = ", bs.list$len))
  print("OLS:")
  print(c(sd.a.ols, sd.b.ols))
  print("GLS:")
  print(c(sd.a.gls, sd.b.gls))
}

bs.25 <- bootstrap.gls(N = 25)
operate.lists(bs.25)
bs.50 <- bootstrap.gls(N = 50)
operate.lists(bs.50)
bs.100 <- bootstrap.gls(N = 100)
operate.lists(bs.100)
bs.500 <- bootstrap.gls(N = 500)
operate.lists(bs.500)
bs.1000 <- bootstrap.gls(N = 1000)
operate.lists(bs.1000)
bs.2000 <- bootstrap.gls(N = 2000)
operate.lists(bs.2000)