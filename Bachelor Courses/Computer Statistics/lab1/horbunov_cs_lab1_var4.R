library(rgl)

set.seed(0)

next.step <- function(i, m, a, c)
{
  # Функція для обчислення наступного члена послідовності.
  # m - модуль
  # a, c - дійсні числа
  # i - попередній член послідовності
  (a*i+c)%%m
}

visual.and.test <- function(I0, m, a, c, N, t = "")
{
  # Створення порожнього масиву з N елементів, 
  # ініціалізація його першого елемента
  I <- vector("numeric", N)
  I[1] <- I0
  
  # Генеруємо наступні N-1 членів послідовності,
  # розмістивши їх у масиві
  for(i in 2:N)
  {
    I[i] <- (a*I[i-1] + c)%%m
  }
  
  # Отримання вектора з дійсними значеннями на [0,1] 
  x <- I/m
  
  # Сортування x (варіаційний ряд)
  xs <- sort(x)
  
  # Побудова графіків емпіричної та теоретичної функції розподілу
  plot(xs, ((1:N)/N), type='s', lwd=2,
       xlim=c(0,1), ylim=c(0,1),
       xlab="t", ylab="F(t)", main=t)
  lines(xs, punif(xs), lty=3, lwd=3, col='red')
  
  # Застосування критерія Колмогорова для перевірки нульововї гіпотези
  # H0 = "P(xi < t) = F_[Unif[0,1]](t)"
  # проти альтернативи
  # H1 = "P(xi < t) != F_[Unif[0,1]](t)"
  print(ks.test(x, "punif", 0, 1, alternative = "greater"))
  
  # Діаграма пар з x
  plot(1:N, x, main=t, cex=0.45, xlab = 'index', ylab = 'value')
  
  # Діаграма трійок з x
  x1 <- x[1:(N-2)]
  x2 <- x[2:(N-1)]
  x3 <- x[3:(N-0)]
  plot3d(x1,x2,x3)
}

park.miller.gen.test <- function()
{
  t.0 <- "Park-Miller Generator:"
  print(t.0)
  visual.and.test(2^15, 2^31-1, 7^5, 0, 500, t.0)
}

custom.generator.test <- function()
{
  t.0 <- "Custom Generator:"
  print(t.0)
  visual.and.test(2^15, 2^31, 75831, 0, 500, t.0)
}

