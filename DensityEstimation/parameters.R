## SIMULATION PARAMETERS -------------------------------------------------------
# Sample Size
N <- 1000

# Number of Processes
processes <- 1

# Smoothing Parameters
lambdas <- 10^seq(from = -3, to = 3, length.out = 10)
#lambdas <- c(0.1, 0.01)

# Side Lengths
L <- c(1, 10) #, 100, 1000)
