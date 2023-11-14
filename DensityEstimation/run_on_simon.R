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
  
  # Folders
  dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/correction"), showWarnings = FALSE)
  dir.create(paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/ricorrection"), showWarnings = FALSE)
  
  # L2-norm of the Error (without and with correction terms)
  err_corr[[l]] <- matrix(nrow = processes, ncol = length(lambdas))
  err_ricorr[[l]] <- matrix(nrow = processes, ncol = length(lambdas))
  
  # Cross-Validation Error
  CVerr_corr[[l]] <- rep(0, length(lambdas))
  CVerr_ricorr[[l]] <- rep(0, length(lambdas))
  
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
      
      FEM_DEPDE_ricorr <- FEM(exp(solution_DEPDE_ricorr$g)/domain_area, FEMbasis)
      solution_ricorr <- eval.FEM(FEM_DEPDE_ricorr, mesh.eval$nodes)
      
      # True Density
      true_density <- dens.func(data[,1], data[,2], xrange, yrange)
      
      # Error
      err_corr[[l]][proc,i] <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis.eval) %*% (solution_corr - true_density)^2)
      err_ricorr[[l]][proc,i] <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis.eval) %*% (solution_ricorr - true_density)^2)
      
      if(proc == 1){
        
        # Plot of the Density Estimate
        plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/correction/correction_lambda_", i)
        plot(FEM(solution_corr, FEMbasis.eval), colormap = "jet.col",
             m = min(true_density), M = max(true_density), filename = plotname)
        
        plotname <- paste0(plotname, "_image")
        image(FEM(solution_corr, FEMbasis.eval), colormap = "jet.col",
              m = min(true_density), M = max(true_density), filename = plotname)
        
        
        # Plot of the Density Estimate
        plotname <- paste0("pictures/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/ricorrection/ricorrection_lambda_", i)
        plot(FEM(solution_ricorr, FEMbasis.eval), colormap = "jet.col",
             m = min(true_density), M = max(true_density), filename = plotname)
        
        plotname <- paste0(plotname, "_image")
        image(FEM(solution_ricorr, FEMbasis.eval), colormap = "jet.col",
              m = min(true_density), M = max(true_density), filename = plotname)
        
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

save(lambdas,  err_corr, err_ricorr, 
     CVerr_corr, CVerr_ricorr, L, 
     file = paste0("output/errors.RData"))
