# Ваше завдання - проаналізувати, з якими чинниками пов'язано відмінність 
# у тривалості життя в різних державах світу, і чи відрізняється цей зв'язок 
# у різних регіонах світу.

# source("qqplotinterval.R")

library(readxl)
library(corrplot)
library(car)
library(gvlma)

dat <- read_xls("world1995.xls")
dat <- data.frame(dat[,-1])
n <- nrow(dat)
d <- ncol(dat)
cnames <- colnames(dat)
cols.choose <- c("POPULATN", "DENSITY", "URBAN", "RELIG_KOD", 
                 "LIFEEXPF", "LIFEEXPM", "LITERACY", "POP_INCR", 
                 "BABYMORT", "GDP_CAP", "REGION_KOD", "CALORIES",
                 "AIDS", "BIRTH_RT", "DEATH_RT", "AIDS_RT", "LOG_GDP",
                 "LG_AIDSR", "B_TO_D", "FERTILTY", "LOG_POP",
                 "CROPGROW", "LIT_MALE", "LIT_FEMA", "CLIMATE_KOD")
for(col.c in cnames)
{
  
  dat[dat[, col.c] == "?", col.c] <- NA
  if(col.c%in%cols.choose)
  {
    dat[, col.c] <- as.numeric(dat[, col.c])
  }
}
nonresponse <- apply(dat[,cols.choose], 2, function(column) {
  sum(is.na(column))
})

corr.pearson <- cor(dat[,cols.choose[-c(4,11,25)]], 
                    method = "pearson", use="complete.obs")
corr.spearman <- cor(dat[,cols.choose[-c(4,11,25)]], 
                     method = "spearman", use="complete.obs")
corrplot(corr.pearson, method = "color", number.cex=0.75,
         type = 'upper', addCoef.col = 'black')
corrplot(corr.spearman, method = "color", number.cex=0.75,
         type = 'upper', addCoef.col = 'black')

par(mfrow = c(2,2))

plot(LIFEEXPM ~ LIT_MALE, data = dat)
plot(LIFEEXPM ~ LOG_GDP, data = dat)
plot(LIFEEXPM ~ BIRTH_RT, data = dat)
plot(LIFEEXPM ~ DEATH_RT, data = dat)

plot(LIFEEXPF ~ LIT_FEMA, data = dat)
plot(LIFEEXPF ~ LOG_GDP, data = dat)
plot(LIFEEXPF ~ BIRTH_RT, data = dat)
plot(LIFEEXPF ~ DEATH_RT, data = dat)

par(mfrow = c(1,1))

# Перша модель, чоловіки

lm.based <- lm(LIFEEXPM ~ I((LIT_MALE)^2) + LOG_GDP + BIRTH_RT + DEATH_RT, data = dat)
print(summary(lm.based))
print(gvlma(lm.based))

fit.lm.based <- fitted(lm.based)
resid.lm.based <- resid(lm.based)

plot(fit.lm.based, resid.lm.based, 
     main = "Прогноз - залишки, модель №1")

qqnorm(resid.lm.based, main = "QQ-діаграма залишків, модель №1")
qqline(resid.lm.based)

plot(lm.based)

w.shap <- shapiro.test(resid.lm.based)
ncv.full <- ncvTest(lm.based)
dw.full <- durbinWatsonTest(lm.based)

# Друга модель, чоловіки

lm.based.red <- lm(LIFEEXPM ~ LOG_GDP + BIRTH_RT + DEATH_RT, data = dat)
print(summary(lm.based.red))
print(gvlma(lm.based.red))

fit.lm.based.red <- fitted(lm.based.red)
resid.lm.based.red <- resid(lm.based.red)

plot(fit.lm.based.red, resid.lm.based.red, 
     main = "Прогноз - залишки, модель №2")

qqnorm(resid.lm.based.red, main = "QQ-діаграма залишків, модель №2")
qqline(resid.lm.based.red)

plot(lm.based.red)

w.shap.red <- shapiro.test(resid.lm.based.red)
ncv.red <- ncvTest(lm.based.red)
dw.red <- durbinWatsonTest(lm.based.red)

# Чи варто робити "стратифікацію" моделі по регіонам?

# OECD = Organization for Economic Co-operation and Development
regions <- unique(dat$REGION)
p <- length(regions)

# Підрахунок залишків в необмеженій моделі
lm.h1.subsets <- list()
for(j in 1:p)
{
  lm.h1.subsets[[j]] <- lm(LIFEEXPM ~ LIT_MALE + LOG_GDP + BIRTH_RT + DEATH_RT,
                           data = subset(dat, REGION == regions[j]))
}

resid.lm.h1 <- unlist(lapply(lm.h1.subsets, function(model) { resid(model) }))

# Підрахунок статистики тесту Фішера

alpha <- 0.05
params.h1 <- 5
F.theor <- qf(1 - alpha, p, n - params.h1 * p)

Sh0 <- sum(resid.lm.based^2)
Sh1 <- sum(resid.lm.h1^2)

F.emp <- ((1 / p) * (Sh0 - Sh1)) / ((1 / (n - 24 - params.h1 * p)) * Sh1)