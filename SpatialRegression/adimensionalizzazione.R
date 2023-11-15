library(fdaPDE)
rm(list=ls())
L = c(1,2,3,4,5)
d = 2
correzione = matrix(0,nrow=length(L), ncol=1)
for(l in 1:length(L)){

N = 10
x = seq(0,L[l], by = L[l]/N)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
mesh = refine.mesh.2D(mesh, maximum_area = L[l]^2/N, minimum_angle = 20)
plot(mesh)
nnodes = dim(mesh$nodes)[1]
print(cat(paste0("nnodes ", nnodes, "\n")))
FEMbasis = create.FEM.basis(mesh)

# Test function
# f = function(x, y, z = 1)
# {
#   coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
#   sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
# }

f = function(x, y)
{
  mvtnorm::dmvnorm(cbind(x,y), mean = c(L[l]/2,L[l]/2), sigma = L[l]^2*diag(2))
}


# Exact solution (pointwise at nodes)
sol_exact = f(mesh$nodes[,1], mesh$nodes[,2])

# Add error to simulate data
set.seed(7893475)
ran = range(sol_exact)
data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
#correzione[l] <- (L[l]^2)^(d-2/d)  # |Omega|^(d-2/d)

correzione[l] <- L[l]^(d-4) 
lambda = 10^seq(-6, 6, length.out=60) / correzione[l] # oppure moltiplica per  Omega^(4/d-1)

options(warn=-1)  
#### Test 6.1: grid with exact GCV
invisible(capture.output(sol<-smooth.FEM(observations=data, 
                                         FEMbasis=FEMbasis, 
                                         lambda=lambda, 
                                         lambda.selection.criterion='grid', 
                                         DOF.evaluation='exact', 
                                         lambda.selection.lossfunction='GCV')))
print(cat(paste0("lambda pos ", sol$optimization$lambda_position, "\n")))
plot(log10(lambda*correzione[l]), sol$optimization$GCV_vector, pch=16, 
     main=paste0("[0,", L[l], "] x [0, ", L[l],"]"))
}
