'%!in%' <- function(x, y) !('%in%'(x, y))

workspace_setup <- function()
{
  workspace <- dirname(sys.frame(1)$ofile)
  
  setwd(workspace)
  
  print(getwd())
}

get_dataframe <- function(path)
{
  result_df <- read.csv(
    path,
    header = F, sep = ','
  )
  
  colnames(result_df) <- c(
    'dat', 'z', 'opn', 'mx', 'mn', 'clo', 'vol'
  )
  
  result_df
}

sample_median <- function(sorted_sample)
{
  result <- median(sorted_sample)
  result
}

iqr_and_quartiles <- function(sample)
{
  sorted_sample <- sort(sample)
  q2 <- sample_median(sorted_sample)
  q1 <- sample_median(sorted_sample[sorted_sample <= q2])
  q3 <- sample_median(sorted_sample[sorted_sample >= q2])
  result <- c(q1, q2, q3, q3 - q1)
  result
}

min_max_and_range <- function(sample)
{
  mx <- max(sample)
  mn <- min(sample)
  result <- c(mn, mx, mx - mn)
}

calc_statistics_for_box_plotting <- function(sample)
{
  result_df <- data.frame(
    t(iqr_and_quartiles(sample)),
    t(min_max_and_range(sample))
  )
  
  names(result_df) <- c(
    "q1", "q2", "q3", "iqr", "mn", "mx", "range"
  )
  
  result_df
}

whiskers <- function(sample)
{
  sorted_sample <- sort(sample)
  
  q2 <- sample_median(sorted_sample)
  
  q1 <- sample_median(sorted_sample[sorted_sample <= q2])
  q3 <- sample_median(sorted_sample[sorted_sample >= q2])
  
  iqr <- q3 - q1
  
  upper_w <- max(sample[sample<=q3+1.5*iqr])
  lower_w <- min(sample[sample>=q1-1.5*iqr])
  
  result <- c(upper_w, lower_w)
  result
}

outliners <- function(sample)
{
  whiskers_values <- whiskers(sample)
  
  outliners_higher <- sample[
    sample > whiskers_values[1]
    ]
  outliners_lower <- sample[
    sample < whiskers_values[2]
    ]
  
  print(
    max(outliners_higher)
    )
  print(
    min(outliners_lower)
    )
  
  result_out <- c(outliners_lower, outliners_higher)
  result_out
}
  
get_lr <- function(data_frame)
{
  lr <- log(data_frame$clo[-nrow(data_frame)]/data_frame$clo[-1])
  lr
}

get_mr <- function(data_frame)
{
  mr <- log(data_frame$mx/data_frame$mn)
  mr
}

get_or <- function(data_frame)
{
  or <- log(data_frame$opn/data_frame$clo)
  or
}

operations <- function(sample)
{
  stats_vals <- calc_statistics_for_box_plotting(sample)
  
  print(stats_vals)
  
  wh <- whiskers(sample)
  
  print('whiskers:')
  
  print(wh)
  
  print('outliners:')
  
  outliners(sample)
}

# ========= STARTUP =========

workspace_setup()

path_to_csv_files <- list.files('./datasets')

for(i in 1:length(path_to_csv_files))
{
  path_to_csv_files[i] = paste('./datasets', path_to_csv_files[i], sep="/")
}

print(path_to_csv_files)

data_frames_list <- list()
titles <- list()

lr_values <- list()
mr_values <- list()
or_values <- list()

for(i in 1:length(path_to_csv_files))
{
  current_path <- path_to_csv_files[i]
  
  split_str <- unlist(strsplit(current_path, "/"))
  
  titles[i] <- split_str[[length(split_str)]]
  
  current_frame <- get_dataframe(current_path)
  
  data_frames_list[[i]] <- current_frame
  
  lr_values[[i]] <- get_lr(data_frames_list[[i]])
  
  mr_values[[i]] <- get_mr(data_frames_list[[i]])
  
  or_values[[i]] <- get_or(data_frames_list[[i]])
}

# ========= CALCULATIONS =========

print('lr, mr, or for the first triplet:')

indices <- c(1, 2, 3)

for(i in indices)
{
  print(titles[i])
  
  print(lr_values[[i]])
  print(mr_values[[i]])
  print(or_values[[i]])
}

# ========= BOXPLOT =========

print('building boxplots:')

# --------- BOXPLOT 1 ---------

for(i in indices)
{
  print(titles[i])
  
  print('lr stats:')
  print(
    operations(
      lr_values[[i]]
    )
  )
  
  print('mr stats:')
  print(
    operations(
      mr_values[[i]]
    )
  )
  
  print('or stats:')
  print(
    operations(
      or_values[[i]]
    )
  )
  
  boxplot(
    lr_values[[i]],
    mr_values[[i]],
    or_values[[i]],
    main = paste(
      "Boxplot for", titles[i],
      sep = " "
    ),
    outline = F
  )
}

# --------- BOXPLOT 2 ---------

bool_value <- as.integer(
  readline(prompt = "Show outliers -> ")
)
  
boxplot(
  main = "lr values",
  lr_values,
  outline = bool_value
)

boxplot(
  main = "mr values",
  mr_values,
  outline = bool_value
)

boxplot(
  main = "or values",
  or_values,
  outline = bool_value
)

# ========= HISTOGRAM =========

print('building histograms:')

# --------- HIST 1 ---------

hist(lr_values[[14]])

hist(mr_values[[5]])

hist(or_values[[17]])

# --------- HIST 2 ---------

num = 16

lr_w_out <- lr_values[[num]]
mr_w_out <- mr_values[[num]]
or_w_out <- or_values[[num]]

lr_w_out <- lr_w_out[lr_w_out%!in%outliners(lr_w_out)]
mr_w_out <- mr_w_out[mr_w_out%!in%outliners(mr_w_out)]
or_w_out <- or_w_out[or_w_out%!in%outliners(or_w_out)]

# ---

intervals <- seq(
  min(
    lr_w_out,
    mr_w_out,
    or_w_out
    ) - 0.01,
  max(
    lr_w_out,
    mr_w_out,
    or_w_out
  ) + 0.01,
  length.out = 75
)

hist(
  main = titles[num],
  lr_w_out, breaks = intervals,
  probability = T,
  col = rgb(1, 0, 0, 0.25),
  ylim=c(0,35)
)

hist(
  mr_w_out, breaks = intervals,
  probability = T,
  col = rgb(0, 1, 0, 0.25),
  add = T
)

hist(
  or_w_out, breaks = intervals,
  probability = T,
  col = rgb(0, 0, 1, 0.25),
  add = T
)

legend(
  x = "topright", y = 20,
  c("lr", "mr", "or"),
  fill=c(
    rgb(1,0,0,0.25),
    rgb(0,1,0,0.25),
    rgb(0,0,1,0.25)
  )
)

# --------- HIST 3 ---------

a_mr <- mr_values[[7]]
b_mr <- mr_values[[13]]
c_mr <- mr_values[[20]]

a_mr <- a_mr[a_mr%!in%outliners(a_mr)]
b_mr <- b_mr[b_mr%!in%outliners(b_mr)]
c_mr <- c_mr[c_mr%!in%outliners(c_mr)]

intervals <- seq(
  min(
    a_mr,
    b_mr,
    c_mr
  ) - 0.01,
  max(
    a_mr,
    b_mr,
    c_mr
  ) + 0.01,
  length.out = 85
)

hist(
  main = "histograms for mr values",
  a_mr, breaks = intervals,
  probability = T,
  col = rgb(1, 0, 0, 0.25),
  ylim=c(0,60)
)

hist(
  b_mr, breaks = intervals,
  probability = T,
  col = rgb(0, 1, 0, 0.25),
  add = T
)

hist(
  c_mr, breaks = intervals,
  probability = T,
  col = rgb(0, 0, 1, 0.25),
  add = T
)

legend(
  x = "topright", y = 20,
  c(titles[7], titles[13], titles[20]),
  fill=c(
    rgb(1,0,0,0.25),
    rgb(0,1,0,0.25),
    rgb(0,0,1,0.25)
  )
)
