workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

c.table <- read.table('c4.txt', header=T, sep='\t', fileEncoding = "UTF16LE")
ln.c.table <- log(c.table)

# Model 1 implementation

lm.log.1 <- lm(Y ~ X, data=ln.c.table)
model.1.coefficients <- lm.log.1$coefficients

# Graphical analysis of Model 1

plot(c.table, cex = 0.5); grid()
curve(exp(model.1.coefficients[1])*x^(model.1.coefficients[2]), col='red', 
      add=T)
#identify(c.table)

plot(ln.c.table, cex = 0.5); grid()
abline(model.1.coefficients, col='red')
#identify(ln.c.table)

m1.residuals <- lm.log.1$residuals
m.r.1 <- mean(m1.residuals[-c(38, 83)])
s.r.1 <- sd(m1.residuals[-c(38, 83)])

# Histogram of residuals, 1

h.1 <- hist(m1.residuals[-c(38, 83)], plot = F) # excluding out-values
delta.1 <- h.1$breaks[2] - h.1$breaks[1]

peak.1 <- dnorm(m.r.1, mean=m.r.1, sd=s.r.1) * length(m1.residuals) * delta.1

plot(h.1, ylim = c(0, max(peak.1, h.1$counts)), density = 15)
curve(length(m1.residuals)*delta.1*dnorm(
  x, mean = m.r.1, sd = s.r.1), 
  col='red', add=T)

# Prediction-residuals, 1

plot(lm.log.1$fitted.values, m1.residuals, cex = 0.5,
     xlab="Predicted", ylab="Residuals"); grid()

# QQ-plot of residuals, 1

qqnorm(m1.residuals, cex=0.5); grid()
qqline(m1.residuals, col='red')

qqnorm(m1.residuals[-c(38, 83)], cex=0.5); grid()
qqline(m1.residuals[-c(38, 83)], col='red')

# Summary of Model 1

print(summary(lm.log.1))

# Model 2 implementation

lm.log.2 <- lm(Y ~ X, data=ln.c.table[-c(38, 83),])
model.2.coefficients <- lm.log.2$coefficients

# Graphical analysis of Model 2

plot(c.table, cex = 0.5); grid()
curve(exp(model.2.coefficients[1])*x^(model.2.coefficients[2]), col='blue', 
      add=T)

plot(ln.c.table, cex = 0.5); grid()
abline(model.2.coefficients, col='blue')

m2.residuals <- lm.log.2$residuals
m.r.2 <- mean(m2.residuals)
s.r.2 <- sd(m2.residuals)

# Histogram of residuals, 2

h.2 <- hist(m2.residuals, plot = F)
delta.2 <- h.2$breaks[2] - h.2$breaks[1]

peak.2 <- dnorm(m.r.2, mean=m.r.2, sd=s.r.2) * length(m2.residuals) * delta.2

plot(h.2, ylim = c(0, max(peak.2, h.2$counts)), density = 15)
curve(length(m2.residuals)*delta.2*dnorm(
  x, mean = m.r.2, sd = s.r.2), 
      col='red', add=T)

# Prediction-residuals, 2

plot(lm.log.2$fitted.values, m2.residuals, cex = 0.5,
     xlab="Predicted", ylab="Residuals"); grid()

# QQ-plot of residuals, 2

qqnorm(m2.residuals, cex=0.5); grid()
qqline(m2.residuals, col='red')

# Summary of Model 2

print(summary(lm.log.2))

