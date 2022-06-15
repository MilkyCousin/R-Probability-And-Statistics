library(car)
library(plotrix)

workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

d <- read.table('kaffee.txt', header=T)[c('dauer', 'alter')]
d$dauer <- log(d$dauer)
b<-boxplot(dauer ~ alter, data = d)
par(mfrow=c(2,3))
for(i in sort(unique(d$alter)))
{
  d.i <- d$dauer[d$alter == i]
  hist(d.i, freq = F, main = paste("dauer <-> alter: level ", i))
  curve(dnorm(x, mean(d.i), sd(d.i)), add = T, col = 'red')
}
par(mfrow=c(1,1))

# mean values - ?

# H0: aj = ai forall i,j
# H1: exists i,: aj != ai

aov.d <- aov(dauer ~ alter, data = d)
print(summary(aov.d))
# 1.66e-10 *** => H1

simult.conf.mean.int <- function(y, f, alpha.0 = 0.05)
{
  M <- length(unique(f))
  alpha <- 1 - (1 - alpha.0)^(1/M)
  n.vect <- tapply(y, f, length)
  m.vect <- tapply(y, f, mean)
  s.vect <- tapply(y, f, sd)
  q.stud <- qt(1 - 0.5*alpha, n.vect - 1)
  i.bound <- q.stud * s.vect/sqrt(n.vect)
  plotCI(
    1:M, y = m.vect, ylab = "Means",
    xaxt = "n", uiw = i.bound
  )
  axis(1, at = 1:M, labels = sort(unique(f)))
}
simult.conf.mean.int(d$dauer, d$alter)
abline(h = mean(d$dauer[d$alter == 1]), col = 'red')
abline(h = (mean(d$dauer[(d$alter == 2)]) + mean(d$dauer[(d$alter == 5)]))/2, 
       col = 'blue')
abline(h = (mean(d$dauer[(d$alter == 3)]) + mean(d$dauer[(d$alter == 4)]))/2, 
       col = 'green')

# variances - ?

# H0: sj = si forall i,j
# H1: exists i,: sj != si

l.x <- leveneTest(d$dauer, as.factor(d$alter))
print(l.x)

# < 2.2e-16 *** => H1

simult.conf.var.int <- function(y, f, alpha.0 = 0.05)
{
  M <- length(unique(f))
  alpha <- 1 - (1 - alpha.0)^(1/M)
  n.vect <- tapply(y, f, length)
  v.vect <- tapply(y, f, var)
  h.1 <- qchisq(0.5 * alpha, df = n.vect - 1)
  h.2 <- qchisq(1 - 0.5 * alpha, df = n.vect - 1)
  plotCI(
    1:M, y = v.vect, ylab = "Variances",
    xaxt = "n", 
    liw = v.vect/h.2 * (n.vect - 1), 
    uiw = v.vect/h.1 * (n.vect - 1)
  )
  axis(1, at = 1:M, labels = sort(unique(f)))
}
simult.conf.var.int(d$dauer, d$alter)