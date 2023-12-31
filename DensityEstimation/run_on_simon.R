library(dplyr)        # for function: between
library(webshot)      # for function: snapshot3d
library(plotly)       # for plots
if(!require(fdaPDE)){
  devtools::install_github(repo="simonepanzeri/fdaPDE") #, ref="master")
}

source("helper_functions.R")
source("plot.R")
source("parameters.R")

# L2-norm of the Error (without and with correction terms)
err_corr <- vector(mode = 'list', length = length(L))
err_ricorr <- vector(mode = 'list', length = length(L))

# Cross-Validation Errors
CVerr_corr <- vector(mode = 'list', length = length(L))
CVerr_ricorr <- vector(mode = 'list', length = length(L))

## SOLUTION WITH CORRECTION TERMS ----------------------------------------------
library("fdaPDE", quietly = TRUE)

print("##### SOLUTION WITH CORRECTION TERMS #####")
for(l in 1:length(L)){
  
  # Range
  xrange <- c(0,L[l])
  yrange <- c(0,L[l])
  
  # L2-norm of the Error (without and with correction terms)
  err_corr[[l]] <- matrix(nrow = processes, ncol = length(lambdas))
  
  # Cross-Validation Error
  CVerr_corr[[l]] <- rep(0, length(lambdas))
  
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
  if(!dir.exists("pictures"))
    dir.create("pictures", showWarnings = FALSE)
  
  if(!dir.exists(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"))){
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]"), showWarnings = FALSE)
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/sample"), showWarnings = FALSE)
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/mesh"), showWarnings = FALSE)
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/true_density"), showWarnings = FALSE)
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/true_density"), showWarnings = FALSE)
  }
  
  if(!dir.exists(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/correction"))){
    dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/correction"), showWarnings = FALSE)
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
      read.table(paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/",N,"data_",proc,".txt"))
      
      # Smoothing Parameter
      lambda <- lambdas[i]
      
      # Solution (with correction terms)
      invisible(capture.output(
        solution_DEPDE_corr <- DE.FEM(data, FEMbasis, lambda, fvec = NULL, heatStep = 0.1, heatIter = 500, 
                                      stepProposals = NULL, tol1 = 1e-5, tol2 = 0, print = FALSE, nfolds = NULL,
                                      nsimulations = 10000, step_method = "Fixed_Step", direction_method = "L-BFGS5",
                                      preprocess_method = "NoCrossValidation", search = "tree")
      ))
      
      invisible(capture.output(
        solution_DEPDE_ricorr <- DE.FEM(data, FEMbasis, lambda*domain_area, fvec = NULL, heatStep = 0.1, heatIter = 500, 
                                        stepProposals = NULL, tol1 = 1e-5, tol2 = 0, print = FALSE, nfolds = NULL,
                                        nsimulations = 10000, step_method = "Fixed_Step", direction_method = "L-BFGS5",
                                        preprocess_method = "NoCrossValidation", search = "tree")
      ))
      
      # Evaluation
      FEM_DEPDE_corr <- FEM(exp(solution_DEPDE_corr$g)/domain_area, FEMbasis)
      solution_corr <- eval.FEM(FEM_DEPDE_corr, mesh.eval$nodes)
      
      # True Density
      true_density <- dens.func(data[,1], data[,2], xrange, yrange)
      
      # Error
      err_corr[[l]][proc,i] <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis.eval) %*% (solution_corr - true_density)^2)
      err_ricorr[[l]][proc,i] <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis.eval) %*% (solution_ricorr - true_density)^2)
      
      if(proc == 1){
        # Plot of the True Estimate
        plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/true_density/true_density")
        if(!file.exists(plotname)){
          plot(FEM(true_density, FEMbasis.eval), colormap = "jet.col",
               m = min(true_density), M = max(true_density), filename = plotname)
        }
        plotname <- paste0(plotname, "_image")
        if(!file.exists(plotname)){
          image(FEM(true_density, FEMbasis.eval), colormap = "jet.col",
                m = min(true_density), M = max(true_density), filename = plotname)
        }
        # Plot of the Sample
        plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/sample/sample")
        if(!file.exists(plotname)){
          plot.sample(coordinates = data, filename = plotname, xrange = xrange, yrange = yrange)
        }
        # Plot of the Mesh
        plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/mesh/mesh")
        if(!file.exists(plotname)){
          plot(mesh, filename = plotname, pch = 19, cex = 0.5)
          plot(mesh.eval, filename = paste0(plotname, ".eval"), pch = 19, cex = 0.5)
        }
        # Plot of the Density Estimate
        plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/correction/correction_lambda_", i)
        if(!file.exists(plotname)){
          plot(FEM(solution_corr, FEMbasis.eval), colormap = "jet.col",
               m = min(true_density), M = max(true_density), filename = plotname)
        }
        plotname <- paste0(plotname, "_image")
        if(!file.exists(plotname)){
          image(FEM(solution_corr, FEMbasis.eval), colormap = "jet.col",
                m = min(true_density), M = max(true_density), filename = plotname)
        }
        t2 <- proc.time()
        
        # Cross-Validation Errors
        invisible(capture.output(
          solution_DEPDE_corr <- DE.FEM(data, FEMbasis, lambdas, fvec = NULL, heatStep = 0.1, heatIter = 500, 
                                        stepProposals = NULL, tol1 = 1e-5, tol2 = 0, print = FALSE, nfolds = 10,
                                        nsimulations = 10000, step_method = "Fixed_Step", direction_method = "L-BFGS5",
                                        preprocess_method = "RightCV", search = "tree")
        ))
        
        CVerr_corr[[l]] <- solution_DEPDE_corr$CV_err
        
        invisible(capture.output(
          solution_DEPDE_ricorr <- DE.FEM(data, FEMbasis, lambdas*domain_area, fvec = NULL, heatStep = 0.1, heatIter = 500, 
                                          stepProposals = NULL, tol1 = 1e-5, tol2 = 0, print = FALSE, nfolds = 10,
                                          nsimulations = 10000, step_method = "Fixed_Step", direction_method = "L-BFGS5",
                                          preprocess_method = "RightCV", search = "tree")
        ))
        
        CVerr_ricorr[[l]] <- solution_DEPDE_ricorr$CV_err
        
        print(paste("CV errors computed in", round((proc.time()-t2)[3],2), "seconds."))
        
      }
      
      print(paste("Process", proc, "for lambda", i, "done in", round((proc.time()-t1)[3],2), "seconds."))
      
    }
    
    print(paste(processes, "processes for lambda", i, "done in", round((proc.time()-t0)[3],2), "seconds."))
    
  }
  
  print(paste("##### L =", L[l], "DONE #####"))
  
}

save(lambdas,  err_corr, 
     CVerr_corr, L, 
     file = paste0("output/errors.RData"))
