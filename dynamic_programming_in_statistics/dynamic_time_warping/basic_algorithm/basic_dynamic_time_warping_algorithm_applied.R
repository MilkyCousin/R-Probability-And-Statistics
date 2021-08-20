workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  setwd(workspace)
  print(getwd())
}

workspace_setup()

if(!exists("basic.dtw", mode = "function")) source("basic_dynamic_time_warping_algorithm_horbunov_iii_cstat.R")

# load data
df.covid <- read.csv('./data_covid/data.csv', sep=',')
colnames(df.covid) <- c('entity', 'code', 'date', 'daily confirmed cases')

# Франція, Іспанія, Британія, Німетчина
idx_fr <- df.covid[,1] == 'France'
idx_sp <- df.covid[,1] == 'Spain'
idx_gb <- df.covid[,1] == 'United Kingdom'
idx_ge <- df.covid[,1] == 'Germany'

df.france <- df.covid[idx_fr,]
df.spain <- df.covid[idx_sp,]
df.britain <- df.covid[idx_gb,]
df.germany <- df.covid[idx_ge,]

s.fr <- df.france[,4]
s.sp <- df.spain[,4]
s.gb <- df.britain[,4]
s.ge <- df.germany[,4]

s.fr.n <- scale(s.fr, center=F)
s.sp.n <- scale(s.sp, center=F)
s.gb.n <- scale(s.gb, center=F)
s.ge.n <- scale(s.ge, center=F)

mxx <- max(c(length(s.fr.n), length(s.sp.n), length(s.gb.n), length(s.ge.n)))

# fr - ge, fr - gb, fr - sp
# ge - gb, ge - sp
# gb - sp

legend.v <- c('France', 'Germany', 'Spain', 'United Kingdom')
colors.v <- c('red', 'orange', 'green', 'blue')

title.n.1 <- 'Daily confirmed COVID-19 cases: before scaling'
title.n.2 <- 'Daily confirmed COVID-19 cases: after scaling'

x.lab <- 'Time'
y.lab <- 'Daily # of cases'

pos.n <- "topleft"

mny.1 <- min(c(s.fr, s.sp, s.gb, s.ge))
mxy.1 <- max(c(s.fr, s.sp, s.gb, s.ge))

plot(s.fr, ylim=c(mny.1, mxy.1), xlim=c(1, mxx), col=colors.v[1], main=title.n.1, type='l', xlab=x.lab, ylab=y.lab)
lines(s.ge, col=colors.v[2])
lines(s.sp, col=colors.v[3])
lines(s.gb, col=colors.v[4])
legend(pos.n, legend=legend.v, col=colors.v, lwd=1)

mny.2 <- min(c(s.fr.n, s.sp.n, s.gb.n, s.ge.n))
mxy.2 <- max(c(s.fr.n, s.sp.n, s.gb.n, s.ge.n))

plot(s.fr.n, ylim=c(mny.2, mxy.2), xlim=c(1, mxx), col=colors.v[1], main=title.n.2, type='l', xlab=x.lab, ylab=y.lab)
lines(s.ge.n, col=colors.v[2])
lines(s.sp.n, col=colors.v[3])
lines(s.gb.n, col=colors.v[4])
legend(pos.n, legend=legend.v, col=colors.v, lwd=1)

fast.dtw.summary <- function(x, y)
{
  res.dtw <- basic.dtw(x, y)
  print(res.dtw$opt.func)
  plot(res.dtw)
}

countries.lst <- list(s.fr.n, s.ge.n, s.sp.n, s.gb.n)
countries.names <- c('France', 'Germany', 'Spain', 'United Kingdom')

for(i in 1:length(countries.lst))
{
  for(j in 1:i)
  {
    if(j != i)
    {
      mny.1 <- min(c(countries.lst[[i]], countries.lst[[j]]))
      mxy.1 <- max(c(countries.lst[[i]], countries.lst[[j]]))
      
      plot(
        countries.lst[[i]], ylim=c(mny.1, mxy.1), xlim=c(1, mxx), col=colors.v[1], 
        main=paste(countries.names[i], countries.names[j], sep=' - '), 
        type='l', xlab=x.lab, ylab=y.lab, lwd=1.5
      )
      
      lines(countries.lst[[j]], col=colors.v[4], lwd=1.5)
      legend(pos.n, legend=c(countries.names[i], countries.names[j]), col=colors.v[c(1,4)], lwd=1)
      
      print(paste('Summary:', paste(countries.names[i], countries.names[j], sep=' - ')))
      fast.dtw.summary(countries.lst[[i]], countries.lst[[j]])
    }
  }
}