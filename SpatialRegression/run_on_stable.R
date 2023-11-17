pacman::p_load("fdaPDE","dplyr", "webshot", "plotly","mvtnorm", "RColorBrewer")
options(warn = -1)
Sys.setenv(OPENSSL_CONF="/dev/null")

source("parameters.R")
#source("../DensityEstimation/helper_functions.R")
source("helper_functions.R")
#source("../DensityEstimation/plot.R")

# L2-norm of the Error (without and with correction terms)
err_nocorr <- vector(mode = 'list', length = length(L))

# Cross-Validation Errors
CVerr_nocorr <- vector(mode= 'list', length=length(L))
CVerr_corr <- vector(mode= 'list', length=length(L))

## SOLUTION WITHOUT CORRECTION TERMS -------------------------------------------
set.seed(0)

cat("##### SOLUTION WITHOUT CORRECTION TERMS #####\n")

# Generate data and plots ------------------------------------------------------
# data generation 
for(l in 1:length(L)){
  # Range
  xrange <- c(0,L[l])
  yrange <- c(0,L[l])
  for(proc in 1:processes){
    generate.data.regr(N, proc, xrange, yrange)
  }
}

# plots
# for(l in 1:length(L)){
#   # Range
#   xrange <- c(0,L[l])
#   yrange <- c(0,L[l])
#   cat("--------- ", "[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"  ," ---------\n")
#   
#   # Mesh for Estimation
#   n <- 20
#   X <- seq(from = 0, to = L[l], by = L[l]/n)
#   Y <- X
#   grid <- expand.grid(X, Y)
#   mesh <- create.mesh.2D(grid)
#   mesh <- refine.mesh.2D(mesh, maximum_area = L[l]^2/n, minimum_angle = 20)
#   FEMbasis <- create.FEM.basis(mesh)
#   
#   # Domain Area
#   domain_area <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis) %*% rep(1,nrow(FEMbasis$mesh$nodes)))
#   
#   # Mesh for Evaluation
#   n <- 25
#   X.eval <- seq(from = 0, to = L[l], by = L[l]/n)
#   Y.eval <- X.eval
#   grid <- expand.grid(X.eval, Y.eval)
#   mesh.eval <- create.mesh.2D(grid)
#   mesh.eval <- refine.mesh.2D(mesh.eval, maximum_area = L[l]^2/n, minimum_angle = 20)
#   FEMbasis.eval <- create.FEM.basis(mesh.eval)
#   
#   # Folders
#   if(!dir.exists("pictures"))
#     dir.create("pictures", showWarnings = FALSE)
#   
#   if(!dir.exists(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"))){
#     dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"), showWarnings = FALSE)
#     dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/sample"), showWarnings = FALSE)
#     dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/mesh"), showWarnings = FALSE)
#     dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/true_density"), showWarnings = FALSE)
#   }
#   
#   if(!dir.exists(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/no_correction"))){
#     dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/no_correction"), showWarnings = FALSE)
#   }else{
#     next
#   }
# 
#     for(i in 1:length(lambdas)){
#       
#       t0 <- proc.time()
#       
#       proc <- 1 
#       # Read the Data
#       data.regr <- read.table(paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/",N,"data_",proc,".txt"))
#       
#       data <- data.regr[,1:2]
#       observations <- data.regr[,4]
#       
#       # Smoothing Parameter
#       lambda <- lambdas[i]
#       
#       #Solution (without correction terms)
#       invisible(capture.output(
#         solution_SRPDE_nocorr <- smooth.FEM(observations=observations, locations = data,
#                                             FEMbasis=FEMbasis,
#                                             lambda=lambda*N)
#       ))
#       
#       # Evaluation
#       FEM_SRPDE_nocorr <- solution_SRPDE_nocorr$fit.FEM
#       solution_nocorr <- eval.FEM(FEM_SRPDE_nocorr, mesh.eval$nodes)
#       
#       # True Density
#       true_density <- dens.func(mesh.eval$nodes , xrange, yrange)
#       
#       # Plot of the True Estimate
#       plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/true_density/true_density.png")
#       if(!file.exists(plotname)){
#         fig <- plotly.density(f = FEM(true_density, FEMbasis.eval),
#                        m = min(true_density), M = max(true_density),
#                        filename = NULL, colorscale = "Jet")
#         save_image(fig, file = plotname)
#       }
#       plotname <- paste0(plotname, "_image")
#       if(!file.exists(plotname)){
#         imagely.density(f = FEM(true_density, FEMbasis.eval),
#                         m = min(true_density), M = max(true_density),
#                         filename = plotname, colorscale = "Jet")
#       }
#       # Plot of the Sample
#       plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/sample/sample")
#       if(!file.exists(plotname)){
#         plotly.sample(coordinates = data, filename = plotname, xrange = xrange, yrange = yrange)
#       }
#       # Plot of the Mesh
#       plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/mesh/mesh")
#       if(!file.exists(plotname)){
#         plotly.mesh(x = mesh, filename = plotname)
#         plotly.mesh(x = mesh.eval, filename = paste0(plotname, ".eval"))
#       }
#       # Plot of the Density Estimate
#       plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/no_correction/no_correction_lambda_", i)
#       if(!file.exists(plotname)){
#         plotly.density(f = FEM(solution_nocorr, FEMbasis.eval),
#                        m = min(true_density), M = max(true_density),
#                        filename = plotname, colorscale = "Jet")
#       }
#       plotname <- paste0(plotname, "_image")
#       if(!file.exists(plotname)){
#         imagely.density(f = FEM(solution_nocorr, FEMbasis.eval),
#                         m = min(true_density), M = max(true_density),
#                         filename = plotname, colorscale = "Jet")
#       }
#     }
# }

# L2-norm of the Error (without correction terms)
err_nocorr <- matrix(nrow = processes, ncol = length(L))
err_corr <- matrix(nrow = processes, ncol = length(L))

pdf("CV_errors.pdf")
for(l in 1:length(L)){
  # Range
  xrange <- c(0,L[l])
  yrange <- c(0,L[l])
  cat("--------- ", "[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"  ," ---------\n")
  
  # Cross-Validation Error
  CVerr_nocorr[[l]] <- matrix(nrow = processes, ncol = length(lambdas))
  CVerr_corr[[l]] <- matrix(nrow = processes, ncol = length(lambdas))
  
  # Mesh for Estimation
  n <- 30
  X <- seq(from = 0, to = L[l], by = L[l]/n)
  Y <- X
  grid <- expand.grid(X, Y)
  mesh <- create.mesh.2D(grid)
  mesh <- refine.mesh.2D(mesh, maximum_area = L[l]^2/n, minimum_angle = 20)
  FEMbasis <- create.FEM.basis(mesh)
  
  # Domain Area
  domain_area <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis) %*% rep(1,nrow(FEMbasis$mesh$nodes)))
  
  # Mesh for Evaluation
  n <- 25
  X.eval <- seq(from = 0, to = L[l], by = L[l]/n)
  Y.eval <- X.eval
  grid <- expand.grid(X.eval, Y.eval)
  mesh.eval <- create.mesh.2D(grid)
  mesh.eval <- refine.mesh.2D(mesh.eval, maximum_area = L[l]^2/n, minimum_angle = 20)
  FEMbasis.eval <- create.FEM.basis(mesh.eval)
  
  for(proc in 1:processes){
    t1 <- proc.time()
    # Generate the Data
    generate.data.regr(N, proc, xrange, yrange)
      
    # Read the Data
    data.regr <- read.table(paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/",N,"data_",
                                   proc,".txt"))
      
    data <- data.regr[,1:2]
    observations <- data.regr[,4]
      
    #Solution (without correction terms)
    invisible(capture.output(
      solution_SRPDE_nocorr <- smooth.FEM(observations=observations, locations=data,
                                          FEMbasis=FEMbasis,
                                          lambda=lambdas, #/domain_area,  #/ domain_area,
                                          lambda.selection.criterion='grid',
                                          DOF.evaluation='exact',
                                          lambda.selection.lossfunction='GCV')
    ))
    
    invisible(capture.output(
      solution_SRPDE_corr <- smooth.FEM(observations=observations, locations=data,
                                          FEMbasis=FEMbasis,
                                          lambda=lambdas/domain_area,  #/ domain_area,
                                          lambda.selection.criterion='grid',
                                          DOF.evaluation='exact',
                                          lambda.selection.lossfunction='GCV')
    ))

    # Evaluation
    FEM_SRPDE_nocorr <- solution_SRPDE_nocorr$fit.FEM
    solution_nocorr <- eval.FEM(FEM_SRPDE_nocorr, mesh.eval$nodes)
    
    FEM_SRPDE_corr <- solution_SRPDE_nocorr$fit.FEM
    solution_corr <- eval.FEM(FEM_SRPDE_corr, mesh.eval$nodes)
    
    # True Density
    true_density <- dens.func(mesh.eval$nodes, xrange, yrange)
    
    # Error
    err_nocorr[proc,l] <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis.eval) %*% (solution_nocorr - true_density)^2)
    err_corr[proc,l] <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis.eval) %*% (solution_corr - true_density)^2)
    CVerr_nocorr[[l]][proc,] <- solution_SRPDE_nocorr$optimization$GCV_vector
    
    CVerr_corr[[l]][proc,] <- solution_SRPDE_corr$optimization$GCV_vector
    plot(log10(lambdas),
         CVerr_nocorr[[l]][proc,], pch=16, type = "l", xlab = "param", ylab = "GCV", lwd=3,
         ylim = c(min(CVerr_nocorr[[l]][proc,], CVerr_corr[[l]][proc,]), 
                  max(CVerr_nocorr[[l]][proc,], CVerr_corr[[l]][proc,])),
         main=paste0("[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"))
    points(log10(lambdas),
         CVerr_corr[[l]][proc,], pch=16, type = "l", xlab = "param", ylab = "GCV", lwd=3, col="red",
         main=paste0("[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"))
    legend("topright", legend = c("no corr", "corr"), lwd=3, col = c("black", "red"))
    

    cat(paste("Process", proc, " done in ", round((proc.time()-t1)[3],2), "seconds.\n"))
  }
}

dev.off()
if(!dir.exists("output"))
  dir.create("output", showWarnings = FALSE)

save(lambdas, err_nocorr,
     CVerr_nocorr, L,
     L, file = paste0("output/errors_nocorr.RData"))

plot(CVerr_nocorr[[l]][proc,])
lambdas
