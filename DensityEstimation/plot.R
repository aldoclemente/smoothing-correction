## GRAPHICAL PARAMETERS --------------------------------------------------------
zoom = 0.8
userMatrix = rbind(c(  1,  0.0000000,  0.0000000,  0),
                   c(  0,  0.7071068,  0.7071068,  0),
                   c(  0, -0.7071068,  0.7071068,  0),
                   c(  0,  0.0000000,  0.0000000,  1))
windowRect = c(50, 50, 400, 400)

#pp <- par3d(no.readonly = TRUE)

## PLOT A FEM OBJECT -----------------------------------------------------------
plot.FEM <- function(x, colormap = "heat.colors", num_refinements = NULL, m = NULL, M = NULL, filename = NULL, ...)
{
  if(is(x$FEMbasis$mesh, "mesh.2D")){
    if(x$FEMbasis$order == 1)
    {
      R_plot.ORD1.FEM(x, colormap, m, M, filename, ...)
    }else{
      R_plot.ORDN.FEM(x, colormap, num_refinements, ...)
    }
  }else if(is(x$FEMbasis$mesh, "mesh.2.5D")){
    R_plot_manifold(x,...)
  }else if(is(x$FEMbasis$mesh, "mesh.3D")){
    R_plot_volume(x,...)
  }else if(is(x$FEMbasis$mesh, "mesh.1.5D")){
    R_plot_graph(x, ...)
  }
}

## PLOT A 2D mesh OBJECT -------------------------------------------------------
plot.mesh.2D <- function(x, filename = NULL, ...)
{

  if(!(is.null(filename))){
    pdf(paste0(filename, ".pdf"), width = 8, height = 8)
  }
  
  par(mai = c(0.15, 0.15, 0.15, 0.15))
  plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", asp = 1, ...)
  segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
           x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
  segments(x$nodes[x$segments[,1],1], x$nodes[x$segments[,1],2],
           x$nodes[x$segments[,2],1], x$nodes[x$segments[,2],2], ...)
  
  if(!(is.null(filename))){
    dev.off()
  }
  
}

## PLOT A 2.5D mesh OBJECT -----------------------------------------------------
plot.mesh.2.5D <- function(x,...){
  
  nodes <- x$nodes
  edges <- as.vector(t(x$edges))
  
  open3d()
  axes3d()
  pop3d("lights")
  light3d(specular="black")
  
  points3d(nodes[,1], nodes[,2], nodes[,3], col="black", ...)
  segments3d(nodes[edges,1], nodes[edges,2], nodes[edges,3], col="black",...)
  
  aspect3d("iso")
  view3d(0,-45)
  
}

## PLOT A 3D mesh OBJECT -------------------------------------------------------
plot.mesh.3D <- function(x,...){
  
  nodes <- x$nodes
  faces <- as.vector(t(x$faces[x$facesmarkers,]))
  
  aux_mesh <- create.mesh.2.5D(nodes=nodes, triangles=x$faces[x$facesmarkers,], order=1)
  edges <- as.vector(t(aux_mesh$edges))
  
  open3d()
  axes3d()
  pop3d("lights")
  light3d(specular="black")
  
  triangles3d(nodes[faces,1],nodes[faces,2],nodes[faces,3],col="white",...)
  points3d(nodes[,1], nodes[,2], nodes[,3], col="black", ...)
  segments3d(nodes[edges,1], nodes[edges,2], nodes[edges,3], col="black",...)
  
  aspect3d("iso")
  view3d(0,-45)
  
}

## PLOT A 1.5D mesh OBJECT -----------------------------------------------------ù
plot.mesh.1.5D <- function(x, ...)
{
  
  if( x$order == 1 ){
    
    plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
    segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
             x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
  }
  else{
    
    plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
    segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
             x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
    points(x$nodes[x$edges[,3],1], x$nodes[x$edges[,3],2],col="red", ...)
  }
}

R_plot.ORD1.FEM = function(FEM, colormap, m = NULL, M = NULL, filename = NULL, ...)
{
  # PLOT  Plots a FEM object FDOBJ over a rectangular grid defined by
  # vectors X and Y;
  #
  
  if (is.null(m)) {m = min(FEM$coeff)}
  if (is.null(M)) {M = max(FEM$coeff)}
  
  nodes <- FEM$FEMbasis$mesh$nodes
  triangles <- as.vector(t(FEM$FEMbasis$mesh$triangles))
  
  colormap <- match.fun(colormap)
  heat <- colormap(100)
  # How many plots are needed?
  nplots <- ncol(FEM$coeff)
  
  for (i in 1:nplots) {
    
    if (i > 1)
      readline("Press any key for the next plot...")
    
    # open3d()
    # axes3d()
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    
    z <- FEM$coeff[triangles,i]
    triangles3d(nodes[triangles,1], nodes[triangles,2], z,
                color = heat[round(99*(z-min(z))/diff(range(z)))+1],...)
    
    # aspect3d(2,2,1)
    # view3d(0,-45)
    aspect3d(2,2,1.25)
    
    if(!(is.null(filename))){
      snapshot3d(paste0(plotname, ".png"), fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      close3d()
    }
  }
}

R_plot.ORDN.FEM = function(FEM, colormap, num_refinements, ...)
{
  # num_refinements sets the number od division on each triangle edge to be applied for rifenment
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  colormap <- match.fun(colormap)
  heat <- colormap(100)
  
  coeff = FEM$coeff
  
  if(is.null(num_refinements))
  {
    num_refinements = 10
  }
  
  # For the reference triangles we construct a regular mesh
  x = seq(from = 0, to = 1, length.out = num_refinements+1)
  y = seq(from = 0, to = 1, length.out = num_refinements+1)
  points_ref = expand.grid(x,y)
  points_ref = points_ref[points_ref[,1] + points_ref[,2] <= 1,]
  
  
  meshi = create.mesh.2D(nodes = points_ref, order = 1)
  #plot(meshi)
  
  # locations is the matrix with that will contain the coordinate of the points where the function is
  # evaluated (1st and 2nd column) and the columns with the evaluation of the ith fucntion on that point
  
  locations = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$nodes), ncol = 2+ncol(coeff))
  triangles = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$triangles), ncol = 3)
  tot = 0
  
  properties<-R_elementProperties(mesh)
  
  for (i in 1:nrow(mesh$triangles))
  {
    # For each traingle we define a fine mesh as the transofrmation of the one constructed for the reference
    transf<-rbind(cbind(properties$transf_coord$diff1x[i],properties$transf_coord$diff2x[i]),c(properties$transf_coord$diff1y[i],properties$transf_coord$diff2y[i]))
    pointsi = t(transf%*%t(meshi$nodes) + mesh$nodes[mesh$triangles[i,1],])
    #We evaluate the fine mesh OBS: we know the triangle we are working on no need for point location
    z = R_eval_local.FEM(FEM, transf=properties, locations = pointsi, element_index = i)
    
    #We store the results
    locations[((i-1)*nrow(pointsi)+1):(i*nrow(pointsi)),] = cbind(pointsi,z)
    triangles[((i-1)*nrow(meshi$triangles)+1):(i*nrow(meshi$triangles)),] = meshi$triangles+tot
    tot = tot + nrow(meshi$nodes)
  }
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d()
    axes3d()
    pop3d("lights")
    light3d(specular="black")
    z = locations[as.vector(t(triangles)), 2 + isurf]
    triangles3d(x = locations[as.vector(t(triangles)) ,1], y = locations[as.vector(t(triangles)) ,2],
                z = z,
                color = heat[round(99*(z-min(z))/(max(z)-min(z)))+1],...)
    aspect3d(2,2,1)
    view3d(0,-45)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

R_plot_manifold = function(FEM, ...){
  nodes <- FEM$FEMbasis$mesh$nodes
  triangles <- as.vector(t(FEM$FEMbasis$mesh$triangles))
  edges <- as.vector(t(FEM$FEMbasis$mesh$edges))
  coeff <- FEM$coeff
  
  p <- jet.col(n=128,alpha=0.8)
  palette(p)
  ncolor <- length(p)
  
  nplots <- ncol(coeff)
  
  for (i in 1:nplots){
    if (i > 1)
      readline("Press any key for the next plot...")
    
    open3d()
    axes3d()
    pop3d("lights")
    light3d(specular="black")
    
    col <- coeff[triangles,i]
    col <- (ncolor-1)*(col-min(col))/diff(range(col))+1
    col <- p[col]
    
    triangles3d(nodes[triangles,1], nodes[triangles,2],
                nodes[triangles,3], color = col,...)
    segments3d(nodes[edges,1], nodes[edges,2], nodes[edges,3],
               color = "black",...)
    aspect3d("iso")
    view3d(0,-45)
  }
}

R_plot_volume = function(FEM,...){
  
  nodes <- FEM$FEMbasis$mesh$nodes
  faces <- as.vector(t(FEM$FEMbasis$mesh$faces[as.logical(FEM$FEMbasis$mesh$facesmarkers),]))
  # edges <- as.vector(t(FEM$FEMbasis$mesh$edges[as.logical(FEM$FEMbasis$mesh$edgesmarkers),]))
  coeff <- FEM$coeff
  
  p <- jet.col(n=128,alpha=0.8)
  palette(p)
  ncolor <- length(p)
  
  nplots <- ncol(FEM$coeff)
  
  for (i in 1:nplots){
    if (i > 1)
      readline("Press any key for the next plot...")
    
    open3d()
    axes3d()
    pop3d("lights")
    light3d(specular="black")
    
    col <- coeff[faces,i]
    col <- (ncolor-1)*(col-min(col))/diff(range(col))+1
    col <- p[col]
    
    triangles3d(nodes[faces,1], nodes[faces,2],
                nodes[faces,3], color = col,...)
    # segments3d(nodes[edges,1], nodes[edges,2], nodes[edges,3],
    #           color = "black",...)
    aspect3d("iso")
    view3d(0,-45)
  }
}

R_plot_graph = function(FEM, ...){
  
  nodes <- FEM$FEMbasis$mesh$nodes
  if(FEM$FEMbasis$order==1){
    edges <- as.vector(t(FEM$FEMbasis$mesh$edges))
  }else{
    edges <- cbind(FEM$FEMbasis$mesh$edges[,1],
                   FEM$FEMbasis$mesh$edges[,3],
                   FEM$FEMbasis$mesh$edges[,3],
                   FEM$FEMbasis$mesh$edges[,2])
    edges <- as.vector(t(edges))
  }
  coeff <- FEM$coeff
  
  p <- jet.col(n=128,alpha=0.8)
  palette(p)
  ncolor <- length(p)
  nplots <- ncol(coeff)
  
  for (i in 1:nplots){
    if (i > 1)
      readline("Press any key for the next plot...")
    open3d()
    pop3d("lights")
    light3d(specular="black")
    
    col <- coeff[edges,i]
    col <- (ncolor-1)*(col-min(col))/diff(range(col))+1
    col <- p[col]
    
    segments3d(nodes[edges,1], nodes[edges,2],rep(0,dim(nodes)[1]),
               color = col,lwd=2.5,...)
    view3d(0,0,zoom=0.75)  
  }
  
}

## IMAGE PLOT OF A 2D FEM OBJECT -----------------------------------------------
image.FEM = function(x, colormap = "heat.colors", num_refinements = NULL, m = NULL, M = NULL, filename = NULL, ...)
{
  if(!is(x$FEMbasis$mesh, "mesh.2D"))
    stop('This function is implemented only for 2D mesh FEM objects')
  
  if(x$FEMbasis$order == 1)
  {
    R_image.ORD1.FEM(x, colormap, m, M, filename, ...)
  }else{
    R_image.ORDN.FEM(x, num_refinements, ...)
  }
}

R_image.ORD1.FEM = function(FEM, colormap, m = NULL, M = NULL, filename = NULL, ...)
{
  # PLOT  Plots a FEM object FDOBJ over a rectangular grid defined by
  # vectors X and Y;
  #
  
  if (is.null(m)) {m = min(FEM$coeff)}
  if (is.null(M)) {M = max(FEM$coeff)}
  
  nodes <- FEM$FEMbasis$mesh$nodes
  triangles <- as.vector(t(FEM$FEMbasis$mesh$triangles))
  
  colormap <- match.fun(colormap)
  heat <- colormap(100)
  # How many plots are needed?
  nplots <- ncol(FEM$coeff)
  
  for (i in 1:nplots) {
    
    if (i > 1)
      readline("Press any key for the next plot...")
    
    #open3d()
    #axes3d()
    open3d(zoom = 0.7, userMatrix = diag(4), windowRect = c(164, 187, 420, 443))
    
    pop3d("lights")
    light3d(specular="black")
    
    z <- FEM$coeff[triangles,i]
    triangles3d(nodes[triangles,1], nodes[triangles,2], 0,
                color = heat[round(99*(z-min(z))/diff(range(z)))+1])
    
    # aspect3d(2,2,1)
    # view3d(0,0)
    aspect3d(2,2,1.25)
    
    if(!(is.null(filename))){
      snapshot3d(paste0(plotname, ".png"), fmt = "png", width = 750, height = 750, webshot = rgl.useNULL())
      close3d()
    }
  }
}

R_image.ORDN.FEM = function(FEM, num_refinements)
{
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  heat = heat.colors(100)
  
  coeff = FEM$coeff
  
  if(is.null(num_refinements))
  {
    num_refinements = 10
  }
  
  x = seq(from = 0, to = 1, length.out = num_refinements+1)
  y = seq(from = 0, to = 1, length.out = num_refinements+1)
  points_ref = expand.grid(x,y)
  points_ref = points_ref[points_ref[,1] + points_ref[,2] <= 1,]
  
  meshi = create.mesh.2D(nodes = points_ref, order = 1)
  #plot(meshi)
  
  locations = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$nodes), ncol = 3*ncol(coeff))
  triangles = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$triangles), ncol = 3*ncol(coeff))
  tot = 0
  
  properties<-R_elementProperties(mesh)
  
  for (i in 1:nrow(mesh$triangles))
  {
    # For each traingle we define a fine mesh as the transofrmation of the one constructed for the reference
    transf<-rbind(cbind(properties$transf_coord$diff1x[i],properties$transf_coord$diff2x[i]),c(properties$transf_coord$diff1y[i],properties$transf_coord$diff2y[i]))
    pointsi = t(transf%*%t(meshi$nodes) + mesh$nodes[mesh$triangles[i,1],])
    #We evaluate the fine mesh OBS: we know the triangle we are working on no need for point location
    z = R_eval_local.FEM(FEM, transf=properties, locations = pointsi, element_index = i)
    
    #We store the results
    locations[((i-1)*nrow(pointsi)+1):(i*nrow(pointsi)),] = cbind(pointsi,z)
    triangles[((i-1)*nrow(meshi$triangles)+1):(i*nrow(meshi$triangles)),] = meshi$triangles+tot
    tot = tot + nrow(meshi$nodes)
  }
  
  heat = heat.colors(100)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d()
    axes3d()
    pop3d("lights")
    light3d(specular="black")
    z = locations[as.vector(t(triangles)), 2 + isurf];
    triangles3d(x = locations[as.vector(t(triangles)) ,1], y = locations[as.vector(t(triangles)) ,2],
                z=0,
                color = heat[round(99*(z- min(z))/(max(z)-min(z)))+1])
    aspect3d(2,2,1)
    view3d(0,0)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

R_eval_local.FEM = function(FEM, transf, locations, element_index)
{
  N = nrow(locations)
  # Augment Xvec and Yvec by ones for computing barycentric coordinates
  Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
  
  # Get nodes and index
  FEMbasis = FEM$FEMbasis
  mesh = FEMbasis$mesh
  nodes = mesh$nodes
  triangles = mesh$triangles
  coeff = FEM$coeff
  nsurf = dim(coeff)[2]
  
  order = FEMbasis$order
  #nodeindex = params$nodeindex
  detJ = transf$detJ
  
  # 1st, 2nd, 3rd vertices of triangles
  
  v1 = nodes[triangles[element_index,1],]
  v2 = nodes[triangles[element_index,2],]
  v3 = nodes[triangles[element_index,3],]
  
  if(order !=2 && order != 1)
  {
    stop('ORDER is neither 1 or 2.')
  }
  
  # Denominator of change of coordinates chsange matrix
  
  modJ = transf$detJ[element_index]
  ones3 = matrix(1,3,1)
  #modJMat = modJ %*% t(ones3)
  
  M1 = c(v2[1]*v3[2] - v3[1]*v2[2], v2[2] - v3[2], v3[1] - v2[1])/(modJ)
  M2 = c(v3[1]*v1[2] - v1[1]*v3[2], v3[2] - v1[2], v1[1] - v3[1])/(modJ)
  M3 = c(v1[1]*v2[2] - v2[1]*v1[2], v1[2] - v2[2], v2[1] - v1[1])/(modJ)
  
  evalmat = matrix(NA, nrow=N, ncol=nsurf)
  
  for (isurf in 1:nsurf)
  {
    for(i in 1:N)
    {
      baryc1 = (M1*Pgpts[i,]) %*% ones3
      baryc2 = (M2*Pgpts[i,]) %*% ones3
      baryc3 = (M3*Pgpts[i,]) %*% ones3
      
      if(order == 2)
      {
        c1 = coeff[triangles[element_index,1],isurf]
        c2 = coeff[triangles[element_index,2],isurf]
        c3 = coeff[triangles[element_index,3],isurf]
        c4 = coeff[triangles[element_index,6],isurf]
        c5 = coeff[triangles[element_index,4],isurf]
        c6 = coeff[triangles[element_index,5],isurf]
        
        fval =  c1*(2* baryc1^2 - baryc1) +
          c2*(2* baryc2^2 - baryc2) +
          c3*(2* baryc3^2 - baryc3) +
          c4*(4* baryc1 * baryc2) +
          c5*(4* baryc2 * baryc3) +
          c6*(4* baryc3 * baryc1)
        evalmat[i,isurf] = fval
      }else{
        c1 = coeff[triangles[element_index,1],isurf]
        c2 = coeff[triangles[element_index,2],isurf]
        c3 = coeff[triangles[element_index,3],isurf]
        fval = c1*baryc1 + c2*baryc2 + c3*baryc3
        evalmat[i,isurf] = fval
      }
    }
  }
  return(evalmat)
}

R_elementProperties=function(mesh)
{
  nele = dim(mesh$triangles)[[1]]
  nodes = mesh$nodes
  triangles = mesh$triangles
  
  #detJ   = matrix(0,nele,1)      #  vector of determinant of transformations
  #metric = array(0,c(nele,2,2))  #  3-d array of metric matrices
  #transf = array(0,c(nele,2,2))
  
  transf_coord = NULL
  transf_coord$diff1x = nodes[triangles[,2],1] - nodes[triangles[,1],1]
  transf_coord$diff1y = nodes[triangles[,2],2] - nodes[triangles[,1],2]
  transf_coord$diff2x = nodes[triangles[,3],1] - nodes[triangles[,1],1]
  transf_coord$diff2y = nodes[triangles[,3],2] - nodes[triangles[,1],2]
  
  #  Jacobian or double of the area of triangle
  detJ = transf_coord$diff1x*transf_coord$diff2y - transf_coord$diff2x*transf_coord$diff1y
  
  #Too slow, computed only for stiff from diff1x,diff1y,..
  # for (i in 1:nele)
  # {
  #   #transf[i,,] = rbind(cbind(diff1x,diff2x),c(diff1y,diff2y))
  #   #  Compute controvariant transformation matrix OSS: This is (tranf)^(-T)
  #   Ael = matrix(c(diff2y, -diff1y, -diff2x,  diff1x),nrow=2,ncol=2,byrow=T)/detJ[i]
  #
  #   #  Compute metric matrix
  #   metric[i,,] = t(Ael)%*%Ael
  # }
  
  #FEStruct <- list(detJ=detJ, metric=metric, transf=transf)
  FEStruct <- list(detJ=detJ, transf_coord=transf_coord)
  return(FEStruct)
}