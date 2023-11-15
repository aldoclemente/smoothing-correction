#'############################################################################'#
#'##################### SPATIO-TEMPORAL DENSITY ESTIMATION ###################'# 
#'###################### HELPER FUNCTIONS FOR SIMULATION 1 ###################'#
#'############################################################################'#

## MIXTURE WEIGHTS (NON-NORMALIZED PRIORS) -------------------------------------
get.priors <- function(){
  # var1 <- mean(eigen(get.parameters(1)[[2]])$values)
  # var2 <- mean(eigen(get.parameters(2)[[2]])$values)
  # var3 <- mean(eigen(get.parameters(3)[[2]])$values)
  # var4 <- mean(eigen(get.parameters(4)[[2]])$values)
  # priors <- c(var1, var2, var3, var4)
  priors <- c(5/20, 5/20, 5/20, 5/20)
  return (priors)
}

## DISTRIBUTION PARAMETERS -----------------------------------------------------
get.parameters <- function(d, xrange, yrange){
  if(d==1){
    # First Gaussian distribution
    mu <- c(-2, -1.5) # Mean
    sigma <- matrix(c(0.8, -0.2, -0.2, 0.8), 2, 2) # Covariance matrix
  }
  else if(d==2){
    # Second Gaussian distribution
    mu <- c(2, -2) # Mean
    sigma <- matrix(c(1.5, 0, 0, 1.5), 2, 2) # Covariance matrix
  }
  else if(d==3){
    # Third Gaussian distribution
    mu <- c(-2, 1.5) # Mean
    sigma <- matrix(c(0.8, 0, 0, 0.8), 2, 2) # Covariance matrix
  }
  else if(d==4){
    # Fourth Gaussian distribution
    mu <- c(2, 2) # Mean
    sigma <- matrix(c(1, 0.9, 0.9, 1), 2, 2) # Covariance matrix
  }
  return (list(mu, sigma))
}

## DATA GENERATION [MIXTURE] ---------------------------------------------------
generate.data.mixture <- function(N, proc, xrange, yrange){
  if(!dir.exists(file.path(getwd(), "data"))){
    dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
    
  }
  if(!dir.exists(file.path(getwd(), 
                           paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/")))){
    dir.create(file.path(getwd(), 
                         paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/")), showWarnings = FALSE)
  
  set.seed(100+proc)
  
  # Number of components of the mixture
  D <- 4 
  
  # Data
  locations <- c()
  distribution <- c()
  n <- 1
  while (n <= N){
    p <- get.priors()
    d <- sample(c(1,2,3,4), 1, replace=T, prob=p)
    
    parameters <- get.parameters(d, xrange, yrange)
    l <- mvtnorm::rmvnorm(n=1, mean=parameters[[1]], sigma=parameters[[2]], method="svd")
    
    if(between(l[1], xrange[1], xrange[2]) & between(l[2], yrange[1], yrange[2])){
      
      distribution <- c(distribution, d)
      
      locations <- rbind(locations, l)
      
      n <- n+1
    }
  }
  data <- cbind(locations)
  
  write.table(data, paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/",N,"data_",proc,".txt"), row.names=F, col.names=F)
  }
}

## DENSITY FUNCTION [MIXTURE] --------------------------------------------------
dens.func.mixture <- function(locations){
  p <- get.priors()
  p <- p/sum(p)
  dens <- p[1]*mvtnorm::dmvnorm(locations, mean=get.parameters(1)[[1]], sigma=get.parameters(1)[[2]]) +
          p[2]*mvtnorm::dmvnorm(locations, mean=get.parameters(2)[[1]], sigma=get.parameters(2)[[2]]) +
          p[3]*mvtnorm::dmvnorm(locations, mean=get.parameters(3)[[1]], sigma=get.parameters(3)[[2]]) +
          p[4]*mvtnorm::dmvnorm(locations, mean=get.parameters(4)[[1]], sigma=get.parameters(4)[[2]])
  
  return (dens)
}

## DATA GENERATION [1 GAUSSIAN] ------------------------------------------------
generate.data <- function(N, proc, xrange, yrange){

  if(!dir.exists(file.path(getwd(), "data"))){
    dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
    
  }
  if(!dir.exists(file.path(getwd(), 
                           paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/"))))
    dir.create(file.path(getwd(), 
                         paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/")), showWarnings = FALSE)
    
  
  if(!file.exists(paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/",N,"data_",proc,".txt"))){
  set.seed(100+proc)
  
  length <- diff(xrange)
  
  data <- mvtnorm::rmvnorm(N, mean = c(length/2, length/2), sigma = 0.025*length^2*diag(2))
  
  write.table(data, paste0("data/[",xrange[1],",",xrange[2],"]x[",yrange[1],",",yrange[2],"]/",N,"data_",proc,".txt"), row.names=F, col.names=F)
  }
}

## DENSITY FUNCTION [1 GAUSSIAN] -----------------------------------------------
dens.func <- function(locations, xrange, yrange){
  
  length <- diff(xrange)
  
  dens <- mvtnorm::dmvnorm(locations, mean = c(length/2, length/2), sigma = 0.025*length^2*diag(2))

  return (dens)
}


## PLOT MESH USING PLOTLY ------------------------------------------------------
plot.mesh <- function(x, filename, ...){ 
  edges <- x$edges
  p <- plot_ly(width = 1000, height = 1000) %>% layout(scene = list(
    aspectratio = list(
      x = 1,
      y = 1
    )),
    xaxis = list(
      title = '',
      showgrid = F,
      zeroline = F,
      showticklabels = F,
      ticks = ''),
    yaxis = list(
      title = '',
      showgrid = F,
      zeroline = F,
      showticklabels = F,
      ticks = ''),
    margin = list(
      b = 0,
      l = 0,
      r = 14,
      t = 13
    )) %>%
    add_markers(x = x$nodes[,1],
                y = x$nodes[,2], 
                color = I('black'),
                hoverinfo = 'text',
                text = paste('</br><b> Coordinates:', round(x$nodes[,1],2),
                             round(x$nodes[,2],2)), 
                showlegend = T, 
                visible = T) %>%
    add_segments(x = x$nodes[edges[,1],1],
                 y = x$nodes[edges[,1],2],
                 xend = x$nodes[edges[,2],1],
                 yend = x$nodes[edges[,2],2], 
                 color = I('black'),
                 showlegend = F) 
  
  plotly::export(p, file = paste0(filename,'.png')) 
  
  # Alternative to Export Images
  # saveWidget(p, paste0(filename, '.html'))
  # webshot(paste0(filename,'.html'), paste0(filename,'.png'), delay = 2)
}

## PLOT SAMPLE USING PLOTLY ----------------------------------------------------
plot.sample <- function(coordinates, filename, xrange, yrange, ...){
  
  DATA <- data.frame(x = coordinates[,1], y = coordinates[,2])
  
  ay <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 2,
    range = yrange
  )
  
  ax <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 2,
    range = xrange
  )
  
  p <- plot_ly(data = DATA, x = ~x, y = ~y, type = 'scatter', mode = 'markers',
               width = 1000, height = 1000,
               marker = list(size = 5,
                             color = 'black',
                             line = list(color = 'black',
                                         width = 2))) %>%
    layout(scene = list(
      aspectmode = "data",
      aspectratio = list(
        x = 1,
        y = 1
      )),
      xaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F,
        ticks = ''),
      yaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F,
        ticks = ''),
      margin = list(
        b = 0,
        l = 0,
        r = 14,
        t = 13
      ))
  
  p <- p %>% layout(xaxis = ax, yaxis = ay)
  
  plotly::export(p, file = paste0(filename,'.png'))
  
  # Alternative to Export Images
  # saveWidget(p, paste0(filename, '.html'))
  # webshot(paste0(filename,'.html'), paste0(filename,'.png'), delay = 2)
}

## PLOT DENSITY USING PLOTLY ---------------------------------------------------
plot.density <- function(f, m = NULL, M = NULL, filename, ...){
  
  if (is.null(m)) {m = min(f$coeff)}
  if (is.null(M)) {M = max(f$coeff)}
  
  DATA <- data.frame(x = f$FEMbasis$mesh$nodes[,1], 
                     y = f$FEMbasis$mesh$nodes[,2],
                     z = f$coeff,
                     coeff = f$coeff)
  
  I = (f$FEMbasis$mesh$triangles[,1]-1)
  J = (f$FEMbasis$mesh$triangles[,2]-1)
  K = (f$FEMbasis$mesh$triangles[,3]-1)
  
  ay <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 2,
    range = range
  )
  
  ax <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 2,
    range = range
  )
  
  p <- plot_ly(DATA, x = ~x, y = ~y, z = ~z, intensity = ~z, color = ~z, type = "contour",
               i = I, j = J, k = K, width = 1000, height = 1000, showscale = F,
               contours = list(
                 start = m,
                 end = M,
                 size = (M-m)/8
                 #showlabels = T
               ), ...
  ) %>%
    layout(scene = list(
      aspectmode = "data",
      aspectratio = list(
        x = 1,
        y = 1
      )),
      xaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F,
        ticks = ''),
      yaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F,
        ticks = ''),
      margin = list(
        b = 0,
        l = 0,
        r = 14,
        t = 13
      )
    )
  
  p <- p %>% layout(xaxis = ax, yaxis = ay)
  
  plotly::export(p, file = paste0(filename,'.png'))
  
  # Alternative to Export Images
  # saveWidget(p, paste0(filename, '.html'))
  # webshot(paste0(filename,'.html'), paste0(filename,'.png'), delay = 2)
  
}

## PLOT DENSITY WITHOUT USING PLOTLY -------------------------------------------
image.density <- function(f, colormap = "heat.colors", m = NULL, M = NULL, filename, ...){
  
  if (is.null(m)) {m = min(f$coeff)}
  if (is.null(M)) {M = max(f$coeff)}
  
  n <- sqrt(nrow(f$FEMbasis$mesh$nodes))
  
  X <- unique(f$FEMbasis$mesh$nodes[,1])
  Y <- unique(f$FEMbasis$mesh$nodes[,2])
  Z <- matrix(data = f$coeff, nrow = n, ncol = n)
  
  colormap <- match.fun(colormap)
  color <- colormap(100)
  
  pdf(paste0(filename, ".pdf"), family = 'serif', width = 8, height = 8)
  #x11(width = 8, height = 8)
  par(mai = c(0.25,0.25,0.25,0.25))
  image2D(x = X, y = Y, z = Z, colkey = FALSE,
          col = color, xlab = "", ylab = "",
          contour = list(nlevels = 5, drawlabels = FALSE),
          zlim = c(m, M), asp = 1, xaxt = "n", yaxt = "n")
  dev.off()
  
}