###*******************************************************************
### Cluster Proximity plot
### Like CLUSION (Strehl and Gosh, 2002), just better
###
### Author: Michael Hahsler


cluproxplot <- function(x, labels = NULL, method = NULL, args = NULL,
  axes = TRUE, average = TRUE, lines = TRUE, silhouettes = TRUE,
  col = gray(seq(from = 0, to = 1, length = 100)), 
  linesCol = "blue", main = "cluproxplot",  ...) {

  dissimMeasure <- attr(x, "method")
  x <- as.matrix(x)

  ### set everything to NULL first
  k <- NULL
  labels.ordered <- NULL
  clusterDissMatrix <- NULL
  usedMethod <- c("N/A","N/A")
  names(usedMethod) <- c("inter-cluster", "intra-cluster")


  dim <- dim(x)[1]
  if(dim != dim(x)[2]) 
  stop("cluproxplot not implemented for non-symmetric matrices!")

  ### set default seriation
  if(is.null(method)) {
    method[1] <- "Murtagh"
  }	
  if(is.na(method[2]) && !is.null(labels)) method[2] <- "Murtagh" 


  ### no seriation
  if(pmatch(method[1], "No seriation", nomatch = FALSE)) { 
    order <- NULL
    labels <- NULL
    silhouettes <- FALSE
    usedMethod[1] <- "No seriation"

    ### seriate whole matrix if no labels are given
  }else if(is.null(labels)) {
    order <- seriation(x, method = method[1], args = args)  
    usedMethod[1] <- attr(order, "method")
    silhouettes <- FALSE
     

    ### seriate clusters for given labels
  }else if (!is.null(labels)){
    if(length(labels) != dim) stop("Number of labels in",
      sQuote("labels"), "does not match dimensions of", sQuote("x"))

    k <- max(labels)
    
    ### seriate with average pairwise dissimilarites between clusters
    clusterDissMatrix <- clusterDissimilarity(x, labels)
    
    if(k>2) {
      clusterOrder <- seriation(clusterDissMatrix, method = method[1], 
	args = args)
      usedMethod[1] <- attr(clusterOrder, "method")
    }else{
     clusterOrder <- 1:k
   }


    
    ### determine order for matrix from cluster order
    order <- c()

    if(pmatch(method[2], "No seriation", nomatch = FALSE)) {
      ### no intra-cluster seriation
      for(i in 1 : k) {
	order <- c(order, which(labels == clusterOrder[i]))
      }
      usedMethod[2] <- "No seriation"
    
    }else{

      ### intra-cluster seriation
    
      if(pmatch(method[2], "Silhouette width", nomatch = FALSE)
        || silhouettes == TRUE) {
	### calculate silhouette values for later use
	sil <- silhouette(labels, x)
      }
    
      for(i in 1 : k) {
	take <- which(labels == clusterOrder[i])
	
	### only seriate for >1 elements
	if(length(take) > 1) {

	  if(pmatch(method[2], "Silhouette width", nomatch = FALSE)) {
	   intraOrder <-  order(sil[take, "sil_width"], decreasing = TRUE)
	   attr(intraOrder, "method") <- "Silhouette width"
	   
	  }else{
	    block <- x[take, take, drop = FALSE] 
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
  }

  
  ### we might need unique labels at some point
  labels.unique <-  unique(labels)

  ### reorder matrix
  if(is.null(order)) xReordered <- x
  else xReordered <- x[order, order]

  
  ### color lower triangle panels with gray values for avg. (dis)similarity
  if(average == TRUE && !is.null(clusterDissMatrix)) {
    for(i in 1 : k) {
      for( j in 1 : k) {

	### check empty clusters
	if(is.na(labels.unique[i])) next
	if(is.na(labels.unique[j])) next
	
	### upper panels stay the unchanged
	
	### do lower panels
	if(i > j) { 
	  xReordered[labels == labels.unique[i], 
	    labels == labels.unique[j]] <- clusterDissMatrix[i, j] 
	}
	
	### do diagonal
	if(i == j) {
	  block <- xReordered[labels == labels.unique[i], 
	    labels == labels.unique[j]]

	  block[lower.tri(block, diag = TRUE)] <-  clusterDissMatrix[i, j]

	  xReordered[labels == labels.unique[i],
	    labels == labels.unique[j]] <- block

	}

      }
    }

  }


  ### make a grid
  def.par <- par(no.readonly = TRUE)
  
  if(silhouettes == FALSE) {
    layout(c(1), widths=c(8), heights=c(8), respect = TRUE) 
  }else{
    layout(cbind(1,2), widths=c(8,4), heights=c(8,8), respect = TRUE) 
  }
  
  #c(bottom, left, top, right)
  # default par(mar=c(5, 4, 4, 2)+0.1)
  par(mar=c(5, 6, 3, 2)+0.1)
   
  
  
  ### plot image rotated back counter clockwise by 90 degrees
  image.default(c(1:dim),c(1:dim),
    t(xReordered)[,dim:1], col = col, axes = FALSE, xlab="", ylab="", 
	main = main, ...)


      
  ### prepare for return value
  cluster.description <- NULL


  ### plot cluster borders if we have labels and order
  if(!is.null(labels)) {

    k <- length(labels.unique)

    clusterWidth <- (tabulate(labels)[labels.unique])
    clusterCuts <- cumsum(clusterWidth)
    
    #print(clusterWidth)
    #print(length(labels))

    ### plot cluster labels as axes
    #if(pmatch(axes, "clusters", nomatch = FALSE)) {
    if(axes == TRUE) {
        clusterCenter <- clusterCuts - clusterWidth / 2
      ### top
      mtext(labels.unique, side = 3, line = 0, at = clusterCenter)      
      ### left (plot is 90 degrees rotated!)
      mtext(labels.unique, side = 2, line = 0, at = dim - clusterCenter, las=1) 
    }

    
    ### handle margin size missing!
    #if(pmatch(axes, "labels", nomatch = FALSE)) {
      #  axis(3, at = 1:length(labels), labels = names(labels), las=2, tick = FALSE)
      #axis(4, at = 1:length(labels), labels = rev(names(labels)), las=1, tick = FALSE)
      #}

    
    if(lines == TRUE){
      ### remove last line
      #clusterCuts <- clusterCuts[-length(clusterCuts)]
      ### plot is 90 degrees rotated!
      abline(h = dim - clusterCuts + 0.5, v = clusterCuts + 0.5, col = linesCol)
    }

  
    # plot box after lines
    box(which = "plot" ,col = linesCol)
  
    ### calculate intra-cluster dissimilarity
    dissimilarity <- diag(clusterDissimilarity(xReordered, labels))

    ### generate description
    cluster.description = data.frame(
      position = c(1 : k),
      labels = labels.unique, 
      size = tabulate(labels)[labels.unique],
      avgDissim = dissimilarity[labels.unique])
  }
  
  if(silhouettes == TRUE) {
    #c(bottom, left, top, right)
    # default par(mar=c(5, 4, 4, 2)+0.1)
    par(mar=c(4, 1, 2, 2)+0.1)
  
    ### get and reorder silhouettes
    s <- sil[order,]
    s <- s[,"sil_width"]

    cluster_order <- cluster.description$labels
    tab <- cumsum(tabulate(labels)[cluster_order])

    space <- rep.int(0, times = length(labels))
    space[tab] <- length(labels) / 100
    space<- rev(space)

    barplot(rev(s), horiz = TRUE, col="gray",
      border=0, space=space, name=NULL, xlab="Silhouette width")
    
  }
  
  ### kill grid
  par(def.par)


  ### remove method attibute from order
  attr(order, "method") <- NULL 

  ### return results 
  invisible(list(order = order, method = usedMethod, k = k,
      dissimMeasure =  dissimMeasure,
      description =  cluster.description))
}




###************************************************************************
### inter and intra cluster dissimilarity matrix from 
### a dissimilarity matrix plus labels

clusterDissimilarity <- function(x, labels) {
  x <- as.matrix(x)

  ### kill self-dissimilarity (which is always 0)
  diag(x) <- NA

  k <- max(labels)
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



