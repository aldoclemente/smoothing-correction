pacman::p_load("fdaPDE","dplyr", "webshot", "plotly","mvtnorm", "RColorBrewer")
options(warn = -1)
Sys.setenv(OPENSSL_CONF="/dev/null")
source("helper_functions.R")
source("plot.R")
source("parameters.R")

# L2-norm of the Error (without and with correction terms)
err_nocorr <- vector(mode = 'list', length = length(L))

# Cross-Validation Errors
CVerr_nocorr <- vector(mode = 'list', length = length(L))

## SOLUTION WITHOUT CORRECTION TERMS -------------------------------------------
library("fdaPDE")

cat("##### SOLUTION WITHOUT CORRECTION TERMS #####\n")
for(l in 1:length(L)){
  
  # Range
  xrange <- c(0,L[l])
  yrange <- c(0,L[l])
  cat("--------- ", "[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"  ," ---------\n")
  # L2-norm of the Error (without correction terms)
  err_nocorr[[l]] <- matrix(nrow = processes, ncol = length(lambdas))
  
  # Cross-Validation Error
  CVerr_nocorr[[l]] <- rep(0, length(lambdas))
  
  # Mesh for Estimation
  n <- 20
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

  # True Density
  true_density <- dens.func(mesh.eval$nodes, xrange, yrange)
  
  # Folders
  # Folders
  if(!dir.exists("pictures"))
    dir.create("pictures", showWarnings = FALSE)
  
  if(!dir.exists(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"))){
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"), showWarnings = FALSE)
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/sample"), showWarnings = FALSE)
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/mesh"), showWarnings = FALSE)
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/true_density"), showWarnings = FALSE)
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/no_correction"), showWarnings = FALSE)
  }else{
    next
  }
  
  if(!dir.exists("output")){
    dir.create("output", showWarnings = FALSE)
  }
  
  for(i in 1:length(lambdas)){
    
    t0 <- proc.time()
    for(proc in 1:processes){
      t1 <- proc.time()
      # Generate the Data
      generate.data(N, proc, xrange, yrange)
      
      # Read the Data
      data <- read.table(paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/",N,"data_",proc,".txt"))
      
      # Smoothing Parameter
      lambda <- lambdas[i]
      
      # Solution (without correction terms)
      invisible(capture.output(
        solution_DEPDE_nocorr <- DE.FEM(data, FEMbasis, lambda, fvec = NULL, heatStep = 0.1, heatIter = 500, 
                                        stepProposals = NULL, tol1 = 1e-5, tol2 = 0, print = FALSE, nfolds = NULL,
                                        nsimulations = 10000, step_method = "Fixed_Step", direction_method = "L-BFGS5",
                                        preprocess_method = "NoCrossValidation", search = "tree")
      ))
      
      # Evaluation
      FEM_DEPDE_nocorr <- FEM(exp(solution_DEPDE_nocorr$g), FEMbasis)
      solution_nocorr <- eval.FEM(FEM_DEPDE_nocorr, mesh.eval$nodes)
      
      # Error
      err_nocorr[[l]][proc,i] <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis.eval) %*% (solution_nocorr - true_density)^2)
      
      if(proc == 1){
        # Plot of the True Estimate
        plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/true_density/true_density")
        plot(FEM(true_density, FEMbasis.eval), colormap = "jet.col",
             m = min(true_density), M = max(true_density), filename = plotname)
        
        plotname <- paste0(plotname, "_image")
        image(FEM(true_density, FEMbasis.eval), colormap = "jet.col",
              m = min(true_density), M = max(true_density), filename = plotname)
        
        # Plot of the Sample
        plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/sample/sample")
        plot.sample(coordinates = data, filename = plotname, xrange = xrange, yrange = yrange)
        
        # Plot of the Mesh
        plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/mesh/mesh")
        plot(mesh, filename = plotname, pch = 19, cex = 0.5)
        plot(mesh.eval, filename = paste0(plotname, ".eval"), pch = 19, cex = 0.5)
        
        # Plot of the Density Estimate
        plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/no_correction/no_correction_lambda_", i)
        plot(FEM(solution_nocorr, FEMbasis.eval), colormap = "jet.col",
             m = min(true_density), M = max(true_density), filename = plotname)
        
        plotname <- paste0(plotname, "_image")
        image(FEM(solution_nocorr, FEMbasis.eval), colormap = "jet.col",
              m = min(true_density), M = max(true_density), filename = plotname)
        
        t2 <- proc.time()
        
        # Cross-Validation Errors
        invisible(capture.output(
          solution_DEPDE_nocorr <- DE.FEM(data, FEMbasis, lambdas, fvec = NULL, heatStep = 0.1, heatIter = 500, 
                                          stepProposals = NULL, tol1 = 1e-5, tol2 = 0, print = FALSE, nfolds = 10,
                                          nsimulations = 10000, step_method = "Fixed_Step", direction_method = "L-BFGS5",
                                          preprocess_method = "RightCV", search = "tree")
        ))
        
        CVerr_nocorr[[l]] <- solution_DEPDE_nocorr$CV_err
        
        cat(paste("CV errors computed in", round((proc.time()-t2)[3],2), "seconds.\n"))
        
      }
      
      cat(paste("Process", proc, "for lambda", i, "done in", round((proc.time()-t1)[3],2), "seconds.\n"))
    }
    cat(paste(processes, "processes for lambda", i, "done in", round((proc.time()-t0)[3],2), "seconds.\n"))
  }
}

save(lambdas, err_nocorr,
     CVerr_nocorr, L,
     L, file = paste0("output/errors_nocorr.RData"))
