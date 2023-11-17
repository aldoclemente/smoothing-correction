## SIMULATION PARAMETERS -------------------------------------------------------
# Sample Size
N <- 500

# Number of Processes
processes <- 1

# Smoothing Parameters
lambdas <- 10^seq(from = -4, to = 2, length.out = 30)
#lambdas <- c(0.1, 0.01)

# Side Lengths
L <- c(1, 10, 100, 1000)