sample_mean <- function(sample)
{
  result <- mean(sample)
  result
}

geometric_mean <- function(sample)
{
  result <- prod(sample) ^ (1/length(sample))
  result
}

harmonic_mean <- function(sample)
{
  result <- length(sample)/sum(1/sample)
  result
}

sample_median <- function(sorted_sample)
{
  result <- median(sorted_sample)
  result
}

mid_range <- function(sample)
{
  result <- 0.5 * (max(sample) + min(sample))
}

variance <- function(sample)
{
  result <- var(sample)
  result
}

standard_deviation <- function(sample)
{
  result <- sd(sample)
  result
}

iq <- function(sorted_sample)
{
  q2 <- sample_median(sorted_sample)
  q1 <- sample_median(sorted_sample[sorted_sample <= q2])
  q3 <- sample_median(sorted_sample[sorted_sample >= q2])
  result <- q3 - q1
  result
}

sample_range <- function(sample)
{
  result <- max(sample) - min(sample)
  result
}

cv <- function(sample)
{
  result <- standard_deviation(sample)/sample_mean(sample)
  result
}

examine <- function(sample, filename)
{
  sorted_sample <- sort(sample)

  result_df <- data.frame(
    sample_mean       (sorted_sample), 
    geometric_mean    (sorted_sample),
    harmonic_mean     (sorted_sample), 
    sample_median     (sorted_sample),
    mid_range         (sorted_sample),
    variance          (sorted_sample), 
    standard_deviation(sorted_sample),
    iq                (sorted_sample),
    sample_range      (sorted_sample),
    cv                (sorted_sample)
  )
  
  names(result_df) <- c(
      "mean", "geometric mean", "harmonic mean",
      "median", "mid-range", "variance",
      "standard deviation", "iq", "range", "cv"
  )
  
  write.csv(
    result_df,
    file=filename
  )
}

fitness_data <- read.csv(
    file="/home/fourier-transform/R/r_proj/countries_fitness.csv",
  header=TRUE, sep=";"
)

tomato_data <- read.csv(
  file="/home/fourier-transform/R/r_proj/countries_tomatoes.csv",
  header=TRUE, sep=";"
)

examine(fitness_data[,2], "/home/fourier-transform/R/r_proj/result_frame_1.csv")
examine(tomato_data[,2], "/home/fourier-transform/R/r_proj/result_frame_2.csv")
