## SIMULATION PARAMETERS -------------------------------------------------------
# Sample Size
N <- 1000

# Number of Processes
processes <- 30

# Smoothing Parameters
lambdas <- 10^seq(from = -5, to = 5, length.out = 50)
#lambdas <- c(0.1, 0.01)

# Side Lengths
L <- c(1, 10, 100, 1000)