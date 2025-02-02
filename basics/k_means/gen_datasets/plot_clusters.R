library(dplyr)
library(ggplot2)

plot_file <- function(filename) {
  data <- read.csv(filename, header = TRUE)
  data %>% ggplot(aes(x = x, y = y, color = cluster)) +
  geom_point() +
  coord_fixed() +
  scale_shape_manual(values=c(2, 3))
  base_name <- basename(filename)
  number <- as.integer(sub(".*clusters([0-9]+).csv", "\\1", base_name))
  out_filename <- sprintf("plot%03d.png", number)
  ggsave(out_filename)
}

main <- function() {
  if (length(commandArgs(trailingOnly = TRUE)) == 1) {
    args <- commandArgs(trailingOnly = TRUE)
  } else {
    stop("Wrong number of arguments, expected 1")
  }
  directory_path <- args

  csv_files <- list.files(path = directory_path, 
                          pattern = "*.csv", full.names = TRUE)
  print(csv_files)
  lapply(csv_files, plot_file)
}

main()