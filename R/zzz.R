# Conway-Maxwell-Poisson normalizing constant
COMP <- function(k, Theta) k * log(Theta[1]) - Theta[2] * lfactorial(k)

# Double poisson normalizing constant
dbl_poisson <- function(k, Theta)
  ifelse(k == 0, 0.5 * log(Theta[2]) - Theta[1] * Theta[2],
    0.5 * log(Theta[2]) - (Theta[1] - k) * Theta[2] - k +
    k * (1 - Theta[2]) * log(k) - lfactorial(k) + Theta[2] * k * log(Theta[1]))
