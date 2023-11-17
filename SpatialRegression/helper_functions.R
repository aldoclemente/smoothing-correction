# dens.func <- function(locations, xrange, yrange){
#   
#   length <- diff(xrange)
#   
#   dens <- 10^2*mvtnorm::dmvnorm(locations, mean = c(length/4, 3*length/4), sigma = 0.05*(length/2)^2*diag(2)) + 
#     10^2*mvtnorm::dmvnorm(locations, mean = c(3*length/4, length/4), sigma = 0.05*(length/2)^2*diag(2))
#   
#   return (dens)
# }

dens.func <- function(locations, xrange, yrange){
  
  length <- diff(xrange)
  
  dens <- sin(2*pi*(length-locations[,1])/length)*sin(2*pi*(length-locations[,2])/length)
  return (dens)
}


generate.data.regr <- function(N, proc, xrange, yrange){
    
    if(!dir.exists(file.path(getwd(), "data"))){
      dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
      
    }
    if(!dir.exists(file.path(getwd(), 
                             paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/"))))
      dir.create(file.path(getwd(), 
                           paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/")), showWarnings = FALSE)
    
      set.seed(100+proc)
      
      length <- diff(xrange)
      data <- matrix(nrow=N, ncol=4)
      data[,1:2] <- cbind(runif(N, min=min(xrange), max=max(xrange)),
                          runif(N, min=min(yrange), max=max(yrange)))
      data[,3] <- dens.func(data[,1:2], xrange, yrange)
      
      data[,4] <- data[,3] + rnorm(nrow(data[,1:2]), sd=0.05*diff(range(data[,3])))
      
      write.table(data, paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/",N,"data_",proc,".txt"), 
                  row.names=F, col.names=F)
}


