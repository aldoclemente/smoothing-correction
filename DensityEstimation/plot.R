## PLOT A MESH OBJECT USING PLOTLY ---------------------------------------------
plotly.mesh <- function(x, filename = NULL, ...){ 
  
  xrange <- range(x$nodes[,1])
  yrange <- range(x$nodes[,2])
  
  edges <- x$edges
  
  ay <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 0,
    range = yrange
  )
  
  ax <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 0,
    range = xrange
  )
  
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
  
  p <- p %>% layout(xaxis = ax, yaxis = ay)
  
  if(!is.null(filename)){
    invisible(capture.output(
      plotly::orca(p, file = paste0(filename,'.pdf'))
    ))
    #plotly::export(p, file = paste0(filename,'.png')) 
  } else {
    p
  }
  
}

## PLOT A SAMPLE USING PLOTLY --------------------------------------------------
plotly.sample <- function(coordinates, filename = NULL, xrange, yrange, ...){
  
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
  
  if(!is.null(filename)){
    invisible(capture.output(
      plotly::orca(p, file = paste0(filename,'.pdf'))
    ))
    #plotly::export(p, file = paste0(filename,'.png')) 
  } else {
    p
  }
  
}

## PLOT DENSITY USING PLOTLY ---------------------------------------------------
plotly.density <-  function(f, m = NULL, M = NULL, filename = NULL, ...){
  
  # Available colorscales (specify with the option: colorscale = "...":
  # 1. Grayscale:           "Greys", "Blacks"
  #
  # 2. Sequential Scales:   "Viridis", "Inferno", "Plasma", "Magma", "Cividis"
  #
  # 3. Diverging Scales:    "RdBu", "Picnic", "Portland", "Jet", "Hot", "Blackbody",
  #                         "Earth", "Electric", "YlOrRd", "YlOrBr", "YlGnBu", "YlGn"
  #
  # 4. Categorical Scales:  "Alphabet", "Aptic", "Darkmint", "Picnic", "Portland",
  #                         "Jet", "Viridis", "Inferno", "Plasma", "Magma", "Cividis"
  #
  # For a custom colorscale:
  #   colorscale = list(c(0, "red"), c(0.5, "green"), c(1, "blue"))
  
  if (is.null(m)) {m = min(f$coeff)}
  if (is.null(M)) {M = max(f$coeff)}
  
  plot_data <- data.frame(X = f$FEMbasis$mesh$nodes[,1], 
                          Y = f$FEMbasis$mesh$nodes[,2],
                          Z = f$coeff,
                          coeff = f$coeff)
  I = (f$FEMbasis$mesh$triangles[,1]-1); J = (f$FEMbasis$mesh$triangles[,2]-1); K = (f$FEMbasis$mesh$triangles[,3]-1)
  
  p <- plot_ly(plot_data, type = "mesh3d", x = ~X, y = ~Y,  z = ~Z, 
               i = I, j = J, k = K,
               intensity = ~coeff, color = ~coeff,
               cmin = m,
               cmax = M,
               contours = list(showlabels = FALSE),
               colorbar = list(title = ""),
               showlegend = FALSE,
               ...)
  
  p <- p %>% layout(
    scene = list(
      xaxis = list(title = '',
                   showgrid = F,
                   zeroline = F,
                   showticklabels = F,
                   ticks = '',
                   showline = FALSE),
      yaxis = list(title = '',
                   showgrid = F,
                   zeroline = F,
                   showticklabels = F,
                   ticks = '',
                   showline = FALSE),
      zaxis = list(title = '',
                   showgrid = F,
                   zeroline = F,
                   showticklabels = F,
                   ticks = '',
                   showline = FALSE),
      aspectratio = list(x = 1, y = 1, z = 0.75),
      camera = list(eye = list(x = 1.25, 
                               y = -1.25, 
                               z = 0.5)),
      paper_bgcolor = "rgba(0,0,0,0)"),
    margin = list(l = 0, r = 0, b = 0, t = 1))
  
  p <- p %>% hide_colorbar()
  
  if(!is.null(filename)){
    invisible(capture.output(
      plotly::orca(p, file = paste0(filename,'.pdf'),
                   width = 600, height = 600)
    ))
    #plotly::export(p, file = paste0(filename,'.png'))
  }
}

## IMAGE PLOT OF A 2D FEM OBJECT USING PLOTLY ----------------------------------
imagely.density <- function(f, m = NULL, M = NULL, filename = NULL, ...){
  
  if (is.null(m)) {m = min(f$coeff)}
  if (is.null(M)) {M = max(f$coeff)}
  
  range = range(f$FEMbasis$mesh$nodes[,1])
  
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
  
  if(!is.null(filename)){
    invisible(capture.output(
      plotly::orca(p, file = paste0(filename,'.pdf'))
    ))
    #plotly::export(p, file = paste0(filename,'.png')) 
  } else {
    p
  }
  
}
