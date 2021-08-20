workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

alpha <- 0.05

h.table <- read.table('house01.txt', header=T)

# Побудуйте регресійну модель залежності між віком голови домогосподарства
# (AGE_HEAD) та повними витратами (TOTALEXP). Чи потрібно враховувати при
# побудові моделі наявність ванни або душу (BATH)?

h.d <- h.table[c("AGE_HEAD", "TOTALEXP", "BATH")]

plot(data.matrix(h.d[1]),data.matrix(h.d[2]), 
     col = c('green', 'purple')[data.matrix(h.d[3])], cex = 0.5);grid()

plot(data.matrix(h.d[1]), log(data.matrix(h.d[2])), 
     col = c('green', 'purple')[data.matrix(h.d[3])], cex = 0.5);grid()

plot(data.matrix(h.d[1])[h.d$BATH == 1], 
     log(data.matrix(h.d[2]))[h.d$BATH == 1], 
     col = 'green', cex = 0.5);grid()

plot(data.matrix(h.d[1])[h.d$BATH == 2], 
     log(data.matrix(h.d[2]))[h.d$BATH == 2], 
     col = 'purple', cex = 0.5);grid()


X.u <- as.matrix(
  data.frame(
    b.0.0 = (h.d[3] == 1) * 1,
    b.1.0 = (h.d[3] == 1) * h.d$AGE_HEAD,
    b.0.1 = (h.d[3] == 2) * 1,
    b.1.1 = (h.d[3] == 2) * h.d$AGE_HEAD
  )
)

n <- nrow(X.u)

Y <- log(data.matrix(h.d$TOTALEXP))

A.u <- t(X.u)%*%X.u
print(A.u)
print(det(A.u))

A.u.inv <- solve(A.u)
b.u <- A.u.inv%*%t(X.u)%*%Y
print(b.u)

X.r <- matrix(cbind(1 + numeric(n), h.d$AGE_HEAD), nrow = n, ncol = 2)

A.r <- t(X.r)%*%X.r
print(A.r)
print(det(A.r))

A.r.inv <- solve(A.r)
b.r <- A.r.inv%*%t(X.r)%*%Y
print(b.r)

d <- ncol(X.u)
d.r <- ncol(X.r)

U.u <- Y - X.u%*%b.u
U.r <- Y - X.r%*%b.r

sq.U.u <- sum(U.u^2)
sq.U.r <- sum(U.r^2)
F.chou.stat <- ((1/d.r) * (sq.U.r - sq.U.u))/((1/(n - d)) * sq.U.u)
F.chou.theo <- qf(1 - alpha, d.r, n - d)
print(c(F.chou.stat, F.chou.theo))

X.r.new <- as.matrix(
  data.frame(
    b.0.0 = (h.d[3] == 1) * 1,
    b.0.1 = (h.d[3] == 2) * 1,
    angle = h.d$AGE_HEAD
  )
)

n <- nrow(X.r.new)

A.r.new <- t(X.r.new)%*%X.r.new
print(A.r.new)
print(det(A.r.new))

A.r.new.inv <- solve(A.r.new)
b.r.new <- A.r.new.inv%*%t(X.r.new)%*%Y
print(b.r.new)

U.r.new <- Y - X.r.new%*%b.r.new
sq.U.r.new <- sum(U.r.new^2)
d.new <- 3

F.fis.stat <- ((1/d.new) * (sq.U.r.new - sq.U.u))/((1/(n - d)) * sq.U.u)
F.fis.theo <- qf(1 - alpha, d.new, n - d)
print(c(F.fis.stat, F.fis.theo))
