# Load necessary libraries
library(dplyr)
library(mvtnorm)
# Load ggplot2 library to create visualizations
library(ggplot2)

# Define function to generate Gaussian data given number of samples (n),
# center point (center), covariance matrix (sigma) and label for the class
generate_gaussian_data <- function(n, center, sigma, label) {
  data <- rmvnorm(n, mean <- center, sigma <- sigma)
  # Convert generated data into a dataframe and rename columns to "x" and "y"
  data <- data.frame(data)
  names(data) <- c("x", "y")
  data <- data %>% mutate(class = factor(label))
  return(data)
}

main <- function()
{
  n <- 5000
  center <- c(5, 5)
  sigma <- matrix(c(1, 0, 0, 1), nrow = 2)

  data1 = generate_gaussian_data(n, center, sigma, 1)

  n <- 5000
  center <- c(1, 1)
  sigma <- matrix(c(1, 0, 0, 1), nrow = 2)

  data2 <- generate_gaussian_data(n, center, sigma, 2)

  n <- 5000
  center <- c(9, 9)
  sigma <- matrix(c(1, 0, 0, 1), nrow = 2)

  data3 <- generate_gaussian_data(n, center, sigma, 3)

  n <- 5000
  center <- c(-9, 9)
  sigma <- matrix(c(1, 0, 0, 1), nrow = 2)

  data4 <- generate_gaussian_data(n, center, sigma, 4)

  # Combine both clusters into a single dataset and add dataset identifier
  data <- bind_rows(data1, data2, data3, data4)
  data$dataset <- "1 - Mixture of Gaussians"

  initial_dataset <- data[, c("x", "y")]
  data %>% ggplot(aes(x = x, y = y, color = class)) +
    geom_point() +
    coord_fixed() +
    scale_shape_manual(values = c(2, 3))

  ggsave("plot.png")
  # Write the initial dataset to a csv file named 'dataset.csv' without row names
  write.csv(initial_dataset, file="dataset.csv", row.names = FALSE)
}

main()
