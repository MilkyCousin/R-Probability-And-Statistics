### Ядерні оцінки щільностей

## Ядро Єпанєчнікова та характеристики

epan.kernel <- function(t)
{
  3/4 * (1 - t^2) * (abs(t) < 1)
}

ep.d.sq <- 3/5
ep.D <- 1/5

## Безпосередньо підрахунок ядерної оцінки

dens.estim <- function(x, h, K)
{
  n <- length(x)
  g <- function(t)
  {
    sum(K((t - x) / h)) / (n * h)
  }
  gv <- function(t)
  {
    sapply(t, g)
  }
}

## Правило Сільвермана для h

h.silv.rule <- function(d.sq, D, n, sigma.est)
{
  sigma.est * ((8 * sqrt(pi) * d.sq) / (3 * D^2 * n))^(1/5)
}

h.silv.simple <- function(d.sq, D, x)
{
  h.silv.rule(d.sq, D, length(x), sd(x))
}

h.silv.improved <- function(d.sq, D, x)
{
  h.silv.rule(d.sq, D, length(x), min(sd(x), IQR(x)/1.34))
}

## CV-вибір для h

f.a <- function(t, a)
{
  ((1 + t)^a - (-1)^a) / a
}

epan.kernel.conv <- function(t)
{ 
  ifelse(abs(t) >= 2, 0, {
    z <- -abs(t)
    (1 - z^2) * f.a(z, 1) + 2 * z * (f.a(z, 2) - f.a(z, 4)) + (z^2 - 2) * f.a(z, 3) + f.a(z, 5)
  })
}

CV.h <- function(h, x)
{
  n <- length(x)
  idx <- 1:n
  double.sum <- sum(sapply(idx, function(j) {
    delta <- (x[idx < j] - x[j])/h
    A <- sum(epan.kernel.conv(delta)) * 9/16
    B <- sum(epan.kernel(delta))
    A / n - 2 * B / (n - 1)
  }))
  (ep.d.sq + 2 * double.sum) / (n * h)
}

h.crossvalid <- function(x, h.min, h.max)
{
  optimize(function(h) { CV.h(h, x) }, c(h.min, h.max))$minimum
}

### Дані

wine.data <- read.csv2("wine.csv", header = T)[,c("Alcogol", "Magnesium", "Site")]
N <- nrow(wine.data)

j <- 1
plot(wine.data$Alcogol[wine.data$Site == j], 
     wine.data$Magnesium[wine.data$Site == j], 
     xlab = "Alcogol", ylab = "Magnesium", main = paste("Site:", j))
grid()
j <- 2
plot(wine.data$Alcogol[wine.data$Site == j], 
     wine.data$Magnesium[wine.data$Site == j], 
     xlab = "Alcogol", ylab = "Magnesium", main = paste("Site:", j))
grid()
j <- 3
plot(wine.data$Alcogol[wine.data$Site == j], 
     wine.data$Magnesium[wine.data$Site == j], 
     xlab = "Alcogol", ylab = "Magnesium", main = paste("Site:", j))
grid()

plot(wine.data$Alcogol, wine.data$Magnesium, type = "n", 
     xlab = "Alcogol", ylab = "Magnesium", main = "Scatter plot")
text(wine.data$Alcogol, wine.data$Magnesium, wine.data$Site,
     col = c("red", "blue", "green")[wine.data$Site])
grid()

set.seed(0)

wine.sites <- unique(wine.data$Site)

p <- 0.8
N.train <- floor(p * N)

idx.train <- sample(1:N, N.train)

wine.data.train <- wine.data[idx.train,]
wine.data.test <- wine.data[-idx.train,]

wine.masks <- sapply(wine.sites, function(m) {
  wine.data.train$Site == m
})
N.sites <- apply(wine.masks, 2, sum)

## Перевірка на некорельованість

test.spearman <- cor.test(wine.data$Alcogol, wine.data$Magnesium, 
                          method = "spearman", exact = F)

## Побудова (наївного) емпірично-баєсівського класифікатора

# Оцінюємо апріорний розподіл номера виноградника
pi.j <- sapply(wine.sites, function(j) { mean(wine.data.train$Site == j) })

# Будуємо оцінки щільностей для кожної групи

# Пілотні оцінки для згладжування
h.j.pilot.alc <- sapply(wine.sites, function(j) { 
  h.silv.improved(ep.d.sq, ep.D, 
                  wine.data.train$Alcogol[wine.data.train$Site == j])
  })
h.j.pilot.mag <- sapply(wine.sites, function(j) { 
  h.silv.improved(ep.d.sq, ep.D, 
                  wine.data.train$Magnesium[wine.data.train$Site == j])
})

# CV-оцінки згладжування
h.j.cv.alc <- sapply(wine.sites, function(j) {
  h.crossvalid(wine.data.train$Alcogol[wine.data.train$Site == j],
               0.1 * h.j.pilot.alc[j], 10 * h.j.pilot.alc[j])
})
h.j.cv.mag <- sapply(wine.sites, function(j) {
  h.crossvalid(wine.data.train$Magnesium[wine.data.train$Site == j],
               0.1 * h.j.pilot.mag[j], 10 * h.j.pilot.mag[j])
})

# Ядерні оцінки щільності на основі CV-оцінки згладжування
alc.dens <- sapply(wine.sites, function(j) {
  dens.estim(wine.data.train$Alcogol[wine.data.train$Site == j],
             h.j.cv.alc[j], epan.kernel)
})
mag.dens <- sapply(wine.sites, function(j) {
  dens.estim(wine.data.train$Magnesium[wine.data.train$Site == j],
             h.j.cv.mag[j], epan.kernel)
})

# Будуємо графіки оцінок щільностей

# Для Alcogol
eps <- 1
I.alc <- seq(min(wine.data.train$Alcogol) - eps, max(wine.data.train$Alcogol) + eps, 0.01)

plot(c(I.alc, I.alc, I.alc), as.numeric(sapply(wine.sites, function(j) {
  (alc.dens[[j]])(I.alc)
})), type = "n", xlab = "Alcogol", ylab = "Density",
main = "Графіки оцінок щільностей за Alcogol")
grid()
for(j in wine.sites)
{
  lines(I.alc, (alc.dens[[j]])(I.alc), col = c("red", "blue", "green")[j])
}
legend("topleft", legend = c("Site = 1", "Site = 2", "Site = 3"),
       col = c("red", "blue", "green"), lwd = 1)

# Для Magnesium
eps <- 1
I.mag <- seq(min(wine.data.train$Magnesium) - eps, max(wine.data.train$Magnesium) + eps, 0.01)

plot(c(I.mag, I.mag, I.mag), as.numeric(sapply(wine.sites, function(j) {
  (mag.dens[[j]])(I.mag)
})), type = "n", xlab = "Magnesium", ylab = "Density",
main = "Графіки оцінок щільностей за Magnesium")
grid()
for(j in wine.sites)
{
  lines(I.mag, (mag.dens[[j]])(I.mag), col = c("red", "blue", "green")[j])
}
legend("topleft", legend = c("Site = 1", "Site = 2", "Site = 3"),
       col = c("red", "blue", "green"), lwd = 1)

# Реалізація класифікатора
g.emp.nbayes <- function(x.alc, x.mag)
{
  fj.alc.emp.x <- c((alc.dens[[1]])(x.alc), (alc.dens[[2]])(x.alc), (alc.dens[[3]])(x.alc))
  fj.mag.emp.x <- c((mag.dens[[1]])(x.mag), (mag.dens[[2]])(x.mag), (mag.dens[[3]])(x.mag))
  which.max(pi.j * fj.alc.emp.x * fj.mag.emp.x)
}

# Прогнозування за наївним емп.-байєс. класифікатором
site.hat.emp.nbayes.train <- apply(wine.data.train, 1, function(u) {
  g.emp.nbayes(u[1], u[2])
})
site.hat.emp.nbayes.test <- apply(wine.data.test, 1, function(u) {
  g.emp.nbayes(u[1], u[2])
})

# Частоти помилок класифікації
L.nb.train <- mean(site.hat.emp.nbayes.train != wine.data.train$Site)
L.nb.test <- mean(site.hat.emp.nbayes.test != wine.data.test$Site)

# Таблиці спряженості
l.nb.train <- table(forecast = site.hat.emp.nbayes.train, real = wine.data.train$Site)
l.nb.test <- table(forecast = site.hat.emp.nbayes.test, real = wine.data.test$Site)

## Емпірично-баєсів класифікатор

# Двовимірна щільність

f.est.2d <- function(h1, h2, x1, x2, n, K)
{ function(y1, y2) { sum(K((y1 - x1)/h1) * K((y2 - x2)/h2)) / (n * h1 * h2) } }

# Двовимірна крос-валідація
CV.h.2d <- function(h1, h2, x1, x2)
{
  n <- length(x1); idx <- 1:n
  double.sum <- sum(sapply(idx, function(j) {
    delta1 <- (x1[idx < j] - x1[j])/h1
    delta2 <- (x2[idx < j] - x2[j])/h2
    A <- sum(epan.kernel.conv(delta1) * epan.kernel.conv(delta2)) * (9/16)^2
    B <- sum(epan.kernel(delta1) * epan.kernel(delta2))
    A / n - 2 * B / (n - 1)
  }))
  (ep.d.sq^2 + 2 * double.sum) / (n * h1 * h2)
}
h.crossvalid.2d <- function(x1, x2, h1_0, h2_0, D = 200)
{ 
  vals.h1 <- 0.1 * h1_0 * ((D-1):0)/(D-1) + 10 * h1_0 * (0:(D-1))/(D-1)
  vals.h2 <- 0.1 * h2_0 * ((D-1):0)/(D-1) + 10 * h2_0 * (0:(D-1))/(D-1)
  CV.grid <- outer(vals.h1, vals.h2, function(h1, h2) { 
    apply(cbind(h1, h2), 1, 
          function(u) { CV.h.2d(u[1], u[2], x1, x2) }
          )
    })
  h1.opt <- vals.h1[which.min(apply(CV.grid, 1, min))]
  h2.opt <- vals.h2[which.min(apply(CV.grid, 2, min))]
  c(h1.opt, h2.opt)
}

# Двовимірні оцінки згладжування за CV
#print("2d CV estimation: start")
#hj.1.2d.cv <- h.crossvalid.2d(x1 = wine.data.train$Alcogol[wine.data.train$Site == 1],
#                              x2 = wine.data.train$Magnesium[wine.data.train$Site == 1],
#                              h1 = h.j.cv.alc[1], h2 = h.j.cv.mag[1])
#print("1st")
#hj.2.2d.cv <- h.crossvalid.2d(x1 = wine.data.train$Alcogol[wine.data.train$Site == 2],
#                              x2 = wine.data.train$Magnesium[wine.data.train$Site == 2],
#                              h1 = h.j.cv.alc[2], h2 = h.j.cv.mag[2])
#print("2nd")
#hj.3.2d.cv <- h.crossvalid.2d(x1 = wine.data.train$Alcogol[wine.data.train$Site == 3],
#                              x2 = wine.data.train$Magnesium[wine.data.train$Site == 3],
#                              h1 = h.j.cv.alc[3], h2 = h.j.cv.mag[3])
#print("3rd")
#print("2d CV estimation: finish")

# Двовимірне правило Сільвермана
d <- 2
n1 <- sum(wine.data.train$Site == 1)
n2 <- sum(wine.data.train$Site == 2)
n3 <- sum(wine.data.train$Site == 3)
hj.1.2d.silv <- (4 / (d + 2))^(1/(d + 4)) * n1^(-1/(d + 4)) * c(sd(wine.data.train$Alcogol[wine.data.train$Site == 1]), sd(wine.data.train$Magnesium[wine.data.train$Site == 1]))
hj.2.2d.silv <- (4 / (d + 2))^(1/(d + 4)) * n1^(-1/(d + 4)) * c(sd(wine.data.train$Alcogol[wine.data.train$Site == 2]), sd(wine.data.train$Magnesium[wine.data.train$Site == 2]))
hj.3.2d.silv <- (4 / (d + 2))^(1/(d + 4)) * n1^(-1/(d + 4)) * c(sd(wine.data.train$Alcogol[wine.data.train$Site == 3]), sd(wine.data.train$Magnesium[wine.data.train$Site == 3]))

f1.est <- f.est.2d(h1 = hj.1.2d.silv[1], h2 = hj.1.2d.silv[2], 
                   x1 = wine.data.train$Alcogol[wine.data.train$Site == 1],
                   x2 = wine.data.train$Magnesium[wine.data.train$Site == 1], 
                   n = n1, K = epan.kernel)

f2.est <- f.est.2d(h1 = hj.2.2d.silv[1], h2 = hj.2.2d.silv[2], 
                   x1 = wine.data.train$Alcogol[wine.data.train$Site == 2],
                   x2 = wine.data.train$Magnesium[wine.data.train$Site == 2], 
                   n = n2, K = epan.kernel)

f3.est <- f.est.2d(h1 = hj.3.2d.silv[1], h2 = hj.3.2d.silv[2], 
                   x1 = wine.data.train$Alcogol[wine.data.train$Site == 3],
                   x2 = wine.data.train$Magnesium[wine.data.train$Site == 3], 
                   n = n3, K = epan.kernel)

# Реалізація класифікатора
g.emp.bayes <- function(x.alc, x.mag)
{
  fj.emp.x <- c(f1.est(x.alc, x.mag), f2.est(x.alc, x.mag), f3.est(x.alc, x.mag))
  which.max(pi.j * fj.emp.x)
}

# Прогнозування за емп.-байєс. класифікатором
site.hat.emp.bayes.train <- apply(wine.data.train, 1, function(u) {
  g.emp.bayes(u[1], u[2])
})
site.hat.emp.bayes.test <- apply(wine.data.test, 1, function(u) {
  g.emp.bayes(u[1], u[2])
})

# Частоти помилок класифікації
L.b.train <- mean(site.hat.emp.bayes.train != wine.data.train$Site)
L.b.test <- mean(site.hat.emp.bayes.test != wine.data.test$Site)

# Таблиці спряженості
l.b.train <- table(forecast = site.hat.emp.bayes.train, real = wine.data.train$Site)
l.b.test <- table(forecast = site.hat.emp.bayes.test, real = wine.data.test$Site)

## Побудова проекційного баєсівського класифікатора

direction <- function(beta)
{
  c(cos(beta), sin(beta))
}

projection <- function(x.vect, dir.vect)
{
  as.numeric(crossprod(x.vect, dir.vect))
}

# Побудова класифікатора
g.proj.bayes <- function(x.proj, f1, f2, f3)
{
  f.proj.x <- c(f1(x.proj), f2(x.proj), f3(x.proj))
  which.max(pi.j * f.proj.x)
}

L.proj <- function(beta)
{
  projected.train <- apply(wine.data.train, 1, function(u) {
    projection(u[-3], direction(beta))
  })
  
  h.s.1 <- h.silv.improved(ep.d.sq, ep.D, projected.train[wine.data.train$Site == 1])
  h.cv.1 <- h.crossvalid(projected.train, 0.1 * h.s.1, 10 * h.s.1)
  
  h.s.2 <- h.silv.improved(ep.d.sq, ep.D, projected.train[wine.data.train$Site == 2])
  h.cv.2 <- h.crossvalid(projected.train, 0.1 * h.s.2, 10 * h.s.2)
  
  h.s.3 <- h.silv.improved(ep.d.sq, ep.D, projected.train[wine.data.train$Site == 3])
  h.cv.3 <- h.crossvalid(projected.train, 0.1 * h.s.3, 10 * h.s.3)
  
  f.1 <- dens.estim(projected.train[wine.data.train$Site == 1], h.cv.1, epan.kernel)
  f.2 <- dens.estim(projected.train[wine.data.train$Site == 2], h.cv.2, epan.kernel)
  f.3 <- dens.estim(projected.train[wine.data.train$Site == 3], h.cv.3, epan.kernel)

  hat.proj <- sapply(projected.train, function(u) {
    g.proj.bayes(u, f.1, f.2, f.3)
  })
  
  mean(hat.proj != wine.data.train$Site)
}

beta.res <- optimize(L.proj, lower = 0, upper = 2 * pi)
beta.opt <- beta.res$minimum

# Лахміття

wine.train.projected.opt <- apply(wine.data.train, 1, function(u) {
  projection(u[-3], direction(beta.opt))
})

h.proj.silv.1 <- h.silv.improved(ep.d.sq, ep.D, wine.train.projected.opt[wine.data.train$Site == 1])
h.proj.cv.1 <- h.crossvalid(wine.train.projected.opt, 0.1 * h.proj.silv.1, 10 * h.proj.silv.1)

h.proj.silv.2 <- h.silv.improved(ep.d.sq, ep.D, wine.train.projected.opt[wine.data.train$Site == 2])
h.proj.cv.2 <- h.crossvalid(wine.train.projected.opt, 0.1 * h.proj.silv.2, 10 * h.proj.silv.2)

h.proj.silv.3 <- h.silv.improved(ep.d.sq, ep.D, wine.train.projected.opt[wine.data.train$Site == 3])
h.proj.cv.3 <- h.crossvalid(wine.train.projected.opt, 0.1 * h.proj.silv.3, 10 * h.proj.silv.3)

f.opt.1 <- dens.estim(wine.train.projected.opt[wine.data.train$Site == 1], h.proj.cv.1, epan.kernel)
f.opt.2 <- dens.estim(wine.train.projected.opt[wine.data.train$Site == 2], h.proj.cv.2, epan.kernel)
f.opt.3 <- dens.estim(wine.train.projected.opt[wine.data.train$Site == 3], h.proj.cv.3, epan.kernel)

g.proj.opt <- function(x.alc, x.mag)
{
  x.proj <- projection(c(x.alc, x.mag), direction(beta.opt))
  g.proj.bayes(x.proj, f.opt.1, f.opt.2, f.opt.3)
}

# Прогнозування за проекційним баєсівським класифікатором
site.hat.proj.bayes.train <- apply(wine.data.train, 1, function(u) {
  g.proj.opt(u[1], u[2])
})
site.hat.proj.bayes.test <- apply(wine.data.test, 1, function(u) {
  g.proj.opt(u[1], u[2])
})

# Частоти помилок класифікації
L.pb.train <- mean(site.hat.proj.bayes.train != wine.data.train$Site)
L.pb.test <- mean(site.hat.proj.bayes.test != wine.data.test$Site)

# Таблиці спряженості
l.pb.train <- table(forecast = site.hat.proj.bayes.train, real = wine.data.train$Site)
l.pb.test <- table(forecast = site.hat.proj.bayes.test, real = wine.data.test$Site)

### Графіки

I.alc <- seq(min(wine.data.train$Alcogol) - eps, max(wine.data.train$Alcogol) + eps, 0.02)
I.mag <- seq(min(wine.data.train$Magnesium) - eps, max(wine.data.train$Magnesium) + eps, 0.25)

z.b <- outer(I.alc, I.mag, function(a, m) {
  apply(cbind(a,m), 1, function(u) { g.emp.bayes(u[1], u[2]) })
})
z.nb <- outer(I.alc, I.mag, function(a, m) {
  apply(cbind(a,m), 1, function(u) { g.emp.nbayes(u[1], u[2]) })
})
z.pb <- outer(I.alc, I.mag, function(a, m) {
  apply(cbind(a,m), 1, function(u) { g.proj.opt(u[1], u[2]) })
})

cols = colorRampPalette(c('#ffcfd4','#cfcfff','#cfcfff','#e6ffea'))(24)

filled.contour(I.alc, I.mag, z.b, col = cols,
               plot.axes = {
                 grid(col = "gray48")
                 text(wine.data$Alcogol, wine.data$Magnesium, wine.data$Site, 
                      col = c("red", "blue", "green")[wine.data$Site])
                 axis(1); axis(2); box()
               }, main = "Емпірично-баєсів класифікатор")


filled.contour(I.alc, I.mag, z.nb, col = cols,
               plot.axes = {
                 grid(col = "gray48")
                 text(wine.data$Alcogol, wine.data$Magnesium, wine.data$Site, 
                      col = c("red", "blue", "green")[wine.data$Site])
                 axis(1); axis(2); box()
               }, main = "Наївний баєсів класифікатор")

filled.contour(I.alc, I.mag, z.pb, col = cols,
               plot.axes = {
                 grid(col = "gray48")
                 text(wine.data$Alcogol, wine.data$Magnesium, wine.data$Site, 
                      col = c("red", "blue", "green")[wine.data$Site])
                 axis(1); axis(2); box()
               }, main = "Проекційний баєсів класифікатор")

### Проекція

plot(wine.data$Alcogol, wine.data$Magnesium, type = "n", asp = 1,
     xlim = c(-40, 20), ylim = c(40, 170),
     xlab = "Alcogol", ylab = "Magnesium", main = "Scatter plot with projections")
text(wine.data$Alcogol, wine.data$Magnesium, wine.data$Site,
     col = c("red", "blue", "green")[wine.data$Site])
grid()

dir.opt <- direction(beta.opt)
abline(a = 0, b = dir.opt[2]/dir.opt[1], col = "black")

proj.point <- function(x,y)
{
  p <- projection(c(x,y), dir.opt)
  dir.opt * p
}

set.seed(1)
ind <- sample(1:N, 20)

proj.xy <- apply(cbind(wine.data$Alcogol[ind], wine.data$Magnesium[ind]), 1, function(u) {
  proj.point(u[1], u[2])
})
segments(wine.data$Alcogol[ind], wine.data$Magnesium[ind], proj.xy[1,], proj.xy[2,],
      col = c("red", "blue", "green")[wine.data$Site[ind]])
points(proj.xy[1,], proj.xy[2,],
       col = c("red", "blue", "green")[wine.data$Site[ind]])

# other

I.alc <- seq(-65, 45, 0.25)
I.mag <- seq(55, 165, 0.25)
z.pb <- outer(I.alc, I.mag, function(a, m) {
  apply(cbind(a,m), 1, function(u) { g.proj.opt(u[1], u[2]) })
})

filled.contour(I.alc, I.mag, z.pb, col = cols, asp = 1, 
               plot.axes = {
                 grid(col = "gray48")
                 text(wine.data$Alcogol, wine.data$Magnesium, wine.data$Site, 
                      col = c("red", "blue", "green")[wine.data$Site])
                 abline(a = 0, b = dir.opt[2]/dir.opt[1], col = "black")
                 segments(wine.data$Alcogol[ind], wine.data$Magnesium[ind], proj.xy[1,], proj.xy[2,],
                          col = c("red", "blue", "green")[wine.data$Site[ind]])
                 points(proj.xy[1,], proj.xy[2,],
                        col = c("red", "blue", "green")[wine.data$Site[ind]])
                 xlim = c(min(I.alc), max(I.alc))
                 ylim = c(min(I.mag), max(I.mag))
                 axis(1); axis(2); box()
               }, main = "Проекційний баєсів класифікатор")
