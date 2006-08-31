###*******************************************************************
### Cluster Visualization by Proximity matrix shading
###
### Author: Michael Hahsler


cluproxplot <- function(x, labels = NULL, 
  method = NULL, args = NULL,
  plot = TRUE, plotOptions = NULL, ...) {

  res <- arrange_proximity_matrix(x, labels = labels,
    method = method, args = args)
  if(plot == TRUE) plot(res, plotOptions)

  invisible(res)
}


print.cluProxMatrix <- function(x, ...) {
  d <- dim(x$x)
  k <- if(!is.null(x$k)) x$k else "N/A (1)"

  cat("Rearranged cluster proximity matrix\n")
  cat("\nmatrix dimensions:",d[1],"x",d[2],"\n")

  cat("\nused seriation methods\n")
  print(x$method)
  cat("\nnumber of clusters k:",k, "\n")
  if(!is.null(x$k)) {
    cat("\ncluster description\n")
    print(x$description)
  }
}

plot.cluProxMatrix <- function(x, plotOptions = NULL, ...) {

  ### read plot options
  ### defaults
  options <- list(clusterLabels = TRUE, 
    averages = TRUE, 
    lines = TRUE, 
    silhouettes = TRUE,
    threshold = NULL,
    main = "Cluster proximity plot",
    col = hcl(h = 0, c = 0, l = seq(5, 95, len = 100)), 
    colorkey = TRUE, 
    linesCol = "black",
    newpage = TRUE, 
    pop = TRUE) 

  
  ### now check and add the plot options
  if(!is.null(plotOptions) && length(plotOptions) != 0) {
    m <- pmatch(names(plotOptions), names(options))

    if(any(is.na(m))) stop(paste("Unknown plot option:", 
	names(plotOptions)[is.na(m)], "\n\t"))
 
    for (i in 1:length(m)) {
      options[[m[i]]] <- unlist(plotOptions[i], use.names = FALSE)
    }
  
  } 
  
  ### clear page
  if(options$newpage) grid.newpage()

  m <- x$x
  k <- x$k
  dim <- dim(m)[1]
  labels <- x$labels
  labels.unique <- unique(labels)

  if(is.null(x$sil)) options$silhouettes <- FALSE
  
  ### color lower triangle panels with gray values for avg. (dis)similarity
  if(options$averages == TRUE 
    && !is.null(x$clusterDissMatrix) 
    && !is.null(labels)) {

    for(i in 1 : k) {
      for( j in 1 : k) {

	### check empty clusters
	if(is.na(labels.unique[i])) next
	if(is.na(labels.unique[j])) next

	### upper panels stay the unchanged

	### do lower panels
	if(i > j) { 
	  m[labels == labels.unique[i], 
	    labels == labels.unique[j]] <- x$clusterDissMatrix[i, j] 
	}

	### do diagonal
	if(i == j) {
	  block <- m[labels == labels.unique[i], 
	    labels == labels.unique[j]]

	  block[lower.tri(block, diag = TRUE)] <-  x$clusterDissMatrix[i, j]

	  m[labels == labels.unique[i],
	    labels == labels.unique[j]] <- block

	}
      }
    }
  }

  ### remove entries > threshold
  
  if(options$silhouettes == FALSE) {
    pushViewport(viewport(layout = grid.layout(6, 3,
	  widths = unit(c(3,1,3), c("lines", "null", "lines")),
	  # title, space, image, space, colorkey, space
	  heights = unit(c(3,3,1,1,1,3), 
	    c("lines", "lines", "null", "lines", "lines", "lines")))))

    main_vp <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    image_vp <- viewport(layout.pos.col = 2, layout.pos.row = 3)
    colorkey_vp <- viewport(layout.pos.col = 2, layout.pos.row = 5)

  }else{
    pushViewport(viewport(layout = grid.layout(6, 5,
	  widths = unit(c(3,3,1,1,3), 
	    c("lines", "null", "lines", "null", "lines")),
	  heights = unit(c(3,3,1,1,1,3), 
	    c("lines", "lines", "null", "lines", "lines", "lines")))))

    main_vp <- viewport(layout.pos.col = 2:4, layout.pos.row = 1)
    image_vp <- viewport(layout.pos.col = 2, layout.pos.row = 3)
    barplot_vp <- viewport(layout.pos.col = 4, layout.pos.row = 3)
    colorkey_vp <- viewport(layout.pos.col = 2, layout.pos.row = 5)

  }


  ### main
  pushViewport(main_vp)
  grid.text(options$main, gp = gpar(fontface = "bold", cex = 1.5))
  upViewport(1)

  pushViewport(image_vp)
  grid_simple_image(m, col = options$col, threshold = options$threshold)
  upViewport(1)

  if(options$colorkey == TRUE){
    pushViewport(colorkey_vp)
    grid_colorkey(0, max(m), col = options$col, threshold = options$threshold)
    upViewport(1)
  }

  ### plot cluster borders if we have labels and order
  if(!is.null(labels)) {

    clusterWidth <- (tabulate(labels)[labels.unique])
    clusterCuts <- cumsum(clusterWidth) + 0.5

    if(options$clusterLabels == TRUE) {
      clusterCenter <- clusterCuts - clusterWidth / 2

      seekViewport("image")
      #grid.text(labels.unique, x = clusterCenter, 
	#	y = unit(-1, "lines"), default.unit="native")

      ### above the plot
      grid.text(labels.unique, x = clusterCenter, 
	y = unit(1, "npc") + unit(1, "lines"), default.unit="native")
      ### left of the plot
      grid.text(labels.unique, x = unit(-1, "lines"),
	y = clusterCenter, default.unit="native")
      upViewport(2)
    }

    if(options$lines == TRUE){
      ### remove last line
      clusterCuts <- clusterCuts[-length(clusterCuts)]

      seekViewport("image")
      for(i in 1:k) {

	grid.lines(
	  x = c(0.5, dim + 0.5), 
	  y = clusterCuts[i], 
	  default.unit="native", gp = gpar(col = options$linesCol))

	grid.lines(
	  x = clusterCuts[i], 
	  y = c(0.5, dim + 0.5), 
	  default.unit="native", gp = gpar(col = options$linesCol))

      }

      ### redraw border
      grid.rect(x = 0,5 * dim, y = 0.5 * dim, width = dim - 1, 
	height =  dim - 1, default.units = "native", gp= gpar(col = "black"))

      upViewport(2)

    }
  }

  if(options$silhouettes == TRUE) {

    ### get and reorder silhouettes
    s <- x$sil[,"sil_width"]

    pushViewport(barplot_vp)
    grid_simple_barplot_horiz(s, xlab = "Silhouette width")
    upViewport(1)

  }

  if (options$pop == TRUE) popViewport(1)
  else upViewport(1)

}



arrange_proximity_matrix <- function(x, labels = NULL, method = NULL, 
  args = NULL, ...) {

  dissimMeasure <- attr(x, "method")
  xm <- as.matrix(x)
  if(class(x) != "dist") x <- as.dist(x)
  
  ### set everything to NULL first
  k <- NULL
  sil <- NULL
  labels.ordered <- NULL
  clusterDissMatrix <- NULL
  usedMethod <- c("N/A","N/A")
  names(usedMethod) <- c("inter-cluster", "intra-cluster")


  dim <- dim(xm)[1]
  if(dim != dim(xm)[2]) 
  stop("cluproxplot not implemented for non-symmetric matrices!")

  ### set default seriation
  if(is.null(method)) {
    method[1] <- "Murtagh"
  }	
  if(is.na(method[2]) && !is.null(labels)) method[2] <- "Murtagh" 


  ### no seriation
  if(pmatch(tolower(method[1]), "no seriation", nomatch = FALSE)) { 
    order <- NULL
    labels <- NULL
    usedMethod[1] <- "No seriation"

    ### seriate whole matrix if no labels are given
  }else if(is.null(labels)) {
    order <- seriation(x, method = method[1], args = args)  
    usedMethod[1] <- attr(order, "method")
     

    ### seriate clusters for given labels
  }else if (!is.null(labels)){
    if(length(labels) != dim) stop("Number of labels in",
      sQuote("labels"), "does not match dimensions of", sQuote("x"))

  
    ### get k
    k <- length(unique(labels))
    
    ### seriate with average pairwise dissimilarites between clusters
    clusterDissMatrix <- clusterDissimilarity(xm, labels)
    
    if(k>2) {
      clusterOrder <- seriation(clusterDissMatrix, method = method[1], 
	args = args)
      usedMethod[1] <- attr(clusterOrder, "method")
    }else{
     clusterOrder <- 1:k
   }


    
    ### determine order for matrix from cluster order
    order <- c()

    if(pmatch(tolower(method[2]), "no seriation", nomatch = FALSE)) {
      ### no intra-cluster seriation
      for(i in 1 : k) {
	order <- c(order, which(labels == clusterOrder[i]))
      }
      usedMethod[2] <- "No seriation"
    
    }else{

      ### intra-cluster seriation

      ### calculate silhouette values for later use
      sil <- silhouette(labels, x)

      for(i in 1 : k) {
	take <- which(labels == clusterOrder[i])

	### only seriate for >1 elements
	if(length(take) > 1) {

	  if(pmatch(tolower(method[2]), "silhouette width", nomatch = FALSE)) {
	   intraOrder <-  order(sil[take, "sil_width"], decreasing = TRUE)
	   attr(intraOrder, "method") <- "Silhouette width"
	   
	  }else{
	    block <- xm[take, take, drop = FALSE] 
	    intraOrder <- seriation(block, method = method[2], args = args) 
	  }

	  order <- c(order, take[intraOrder])
	
	}else{
	    order <- c(order, take)
	}

      }
      usedMethod[2] <- attr(intraOrder, "method")
    }


    ### reorder clusterDissMatrix for later
    clusterDissMatrix  <- 
    clusterDissMatrix[clusterOrder, clusterOrder]

    ### prepare order for labels 
    labels <- labels[order]
  
    ### we might need unique labels at some point
    labels.unique <-  unique(labels)
    
  
  }

  ### reorder matrix
  #if(is.null(order)) xReordered <- x
  #else xReordered <- x[order, order]
  if(!is.null(order)) xm <- xm[order, order]
  
  ### prepare for return value
  cluster.description <- NULL

  ### plot cluster borders if we have labels and order
  if(!is.null(labels)) {
    ### reorder silhouettes
    sil <- sil[order,]
    
    ### calculate avg silhouettes
    avgSil <- sapply(labels.unique, function(x) 
      mean(sil[sil[,"cluster"]==x, "sil_width"])) 
    
    ### calculate intra-cluster dissimilarity
    #avgDissim <- diag(clusterDissimilarity(xm, labels))
  
    ### generate description
    cluster.description = data.frame(
      position = c(1 : k),
      label = labels.unique, 
      size = tabulate(labels)[labels.unique],
      #avgDissimilarity = avgDissim[labels.unique],
      avgDissimilarity = diag(clusterDissMatrix)[labels.unique],
      avgSilhouetteWidth = avgSil)
  }
  
  ### remove method attibute from order
  attr(order, "method") <- NULL 

  result <- list(x = xm, 
    labels = labels, 
    method = usedMethod, 
    k = k, 
    clusterDissMatrix =  clusterDissMatrix,
    sil = sil,
    order = order, 
    dissimMeasure =  dissimMeasure,
    description =  cluster.description)
  class(result) <- "cluProxMatrix" 
	      
  invisible(result)
}


###************************************************************************
### inter and intra cluster dissimilarity matrix from 
### a dissimilarity matrix plus labels

clusterDissimilarity <- function(x, labels) {
  if(class(x) != "matrix") x <- as.matrix(x)

  ### kill self-dissimilarity (which is always 0)
  diag(x) <- NA

  k <- length(unique(labels))
  dissMatrix <- matrix(nrow = k, ncol = k)

  ### calculate avg dissimilarity between clusters
  for(i in 1:k) {
    slice <- x[labels == i, , drop = FALSE]
    for(j in 1:i) {
      block <- slice[,labels == j, drop = FALSE]
      val <- mean(as.vector(block), na.rm = TRUE)

      ### fix for clusters of size 1
      if(is.nan(val)) val <- 0

      dissMatrix[i, j] <- val
      dissMatrix[j, i] <- val
    }
  }

  return(dissMatrix)
}



###************************************************************************
### grid helpers

grid_simple_image <- function(x, name = "image", 
  col = gray.colors(12, 0, 1), threshold = NULL) {

  n <-  ncol(x)
  m <-  nrow(x)
  max_x <- max(x)
 
  if(!is.null(threshold)) x[x>threshold] <- NA
  
  div <- 1/length(col)
  
  ### create a viewport
  vp <- viewport(
    xscale = c(0,(n+1)), yscale = c((m+1),0),
    default.unit="native", name = name)
  pushViewport(vp)

  ### make shure we have a color for the maximal value (see floor +1)
  col[length(col)+1] <- col[length(col)]
  
  ### the highest value is lightest color!
  xs <- sapply(c(1:m), "rep.int", times = n)
  grid.rect(x = xs, y = c(1:n), 1, 1, 
    gp = gpar(fill = col[floor(x/max_x/div)+1], col=0), 
    default.units = "native")

  ### make border
  grid.rect(x = (n+1)/2, y = (m+1)/2, width = n, height = m, 
    default.units = "native")

  upViewport(1)
}


grid_simple_barplot_horiz <- function(height, name = "barplot", xlab="") {
  n <-  length(height)

  ### these plots always start at x = 0 or below!
  lim <- c(min(c(height, 0)), max(height))

  ### create a viewport
  vp <- viewport(
    xscale = lim , yscale = c((n+1),0), default.unit="native", name = name)
  pushViewport(vp)

  grid.rect(x = 0, y = 1:n, width = height, height = 1,
    just = c("left", "center"), default.units = "native",
    gp = gpar(col = 0, fill = "lightgray"))

  # hopefuly there is space outside for axes
  grid.xaxis()
  #grid.lines(x = c(0, 0), y = c(0, n+1), default.units = "native")
  grid.text(xlab, y = unit(-3, "lines"))

  upViewport(1)
}

grid_colorkey <- function(min_x, max_x, col, threshold = NULL, 
  name = "colorkey") {
  
  vp <- viewport(
    xscale = c(min_x, max_x), yscale = c(0,1), 
    default.unit="native", name = name)
  pushViewport(vp)

  range <- max_x - min_x
  n <- length(col) 
  width <- range/(n)
  xs <- seq(min_x + width/2, max_x - width/2, length.out = n)
  
  ### dont display the part above the threshold 
  col[xs > threshold] <- NA
  
  grid.rect(x = xs, y = 0, width = width, height = 1,
          just = c("centre", "bottom"), default.units = "native",
    	    gp = gpar(col = 0, fill = col))

  
  ### box
  grid.rect(x = 0, y = 0, width = 1, height = 1,
    just = c("left", "bottom"), default.units = "npc",
    gp = gpar(col = "black"))
  
 grid.xaxis()
  
  upViewport(1)
}
