###*****************************************************************
### Seriation heuristics for symmetric dissimilarity matrices
###
### Author: Michael Hahsler

### calculate order for columns/rows group high values around the diagonal 
### of symmetric matrix x  
seriation <- function(x, method = NULL, args = NULL) {

  ### check built-in methods
  methods <- c(
    "Murtagh (1985), Algorithm B", 
    "Hierarchical Clustering", 
    "Gruvaeus and Wainer (1972)",
    "Optimal Leaf Ordering")

  ### standard seriation is Murtagh
  if(is.null(method)) methodNr <- 1
  else methodNr <- pmatch(method, methods)
  if(is.na(methodNr)) stop (paste("Unknown method:",sQuote(method)))

  ### we need a matrix
  #x <- as.matrix(x) 
  if(class(x) != "dist") x <- as.dist(x)

  
  if(methodNr == 1) {
    order <- .seriationMurtagh(1-x, args)
  }else if (methodNr == 2) {
    order <- .seriationHclust(x, args)
  }else if (methodNr == 3) {
    order <- .seriationGruvaeus(x, args)
  }else if (methodNr == 4) {
    order <- .seriationOptimal(x,args)
  }

  attr(order, "method") <- methods[methodNr]
  return(order)
}

#####################################################################
# seriation methods

# hclust wrapper
.HclustWrapper <- function(x, args) {

    if(class(args$hclust) == "hclust") {
	hclust <- args$hclust  
    } else {
	method <- if(!is.null(args$method)) args$method else "average"
	hclust <- hclust(x, method = method)
    }

    return(hclust)
}


##################################################################
# just hclust vor seriation
.seriationHclust <- function(x, args) {
    return(.HclustWrapper(x, args)$order) 
}

###################################################################
# Gruvaeus, G. and Wainer, H. (1972), "Two Additions to 
# Hierarchical Cluster Analysis", British
# Journal of Mathematical and Statistical Psychology, 25, 200-206.
.seriationGruvaeus <- function(x, args) {

    hclust <- .HclustWrapper(x, args) 

    if(suppressWarnings(require("gclus", quietly = TRUE)) == FALSE)
    stop(paste("Seriation method", sQuote("Gruvaeus and Wainer, 1972"),
	    "requires package", sQuote("gclus."),
	    "Package not installed!"))

    order <- reorder.hclust(hclust, x)$order
    return(order)
}

#########################################################################
# Optimal leave ordering
#  Z. Bar-Joseph, E. D. Demaine, D. K. Gifford, and T. Jaakkola.
#  (2001) Fast Optimal Leaf Ordering for Hierarchical Clustering.
#  Bioinformatics, Vol. 17 Suppl. 1, pp. 22-29.
.seriationOptimal <- function(x, args) {
  
  hclust <- .HclustWrapper(x, args) 

  return(order.optimal(x, hclust$merge)$order) 

}

######################################################################
# Algorithm B
#  F. Murtagh (1985). Multidimensional Cluster Algorithms. Lectures
#  in Computational Statistics, Physica Verlag, pp. 15.
.seriationMurtagh <- function(x, args) {
  ### args not used!

  ### pre-calculate criterion for column pairs in x 
  criterion <- crossprod(as.matrix(x))
  #criterion <- x
  
  return(order.greedy(1-as.dist(criterion))$order)
}
  
### replaced by the C code order.greedy
.seriationMurtagh_R <- function(x, args) {
  ### args not used!
  xm <- as.matrix(x)
  
  k <- dim(xm)[1]

  ### pre-calculate criterion for column pairs in x 
  criterion <- crossprod(xm)

  ### make enough space (0) means no column
  seriation <- integer(k * 2 - 1)
  ### upper and lower already used position
  left.border <- k
  right.border <- k

  ### Step 1: select an arbitrary column and place in the center
  chosen <- sample(1:k, size = 1)
  seriation[k] <- chosen
  placed <- chosen
  criterion[, chosen] <- NA

  ### Step 2: place the column which gives the best criterion
  while(length(placed) < k){
    ### find optimal column
    right.max <- as.integer(which.max(criterion[seriation[right.border],]))
    right.val <- criterion[seriation[right.border],right.max]
    left.max <- as.integer(which.max(criterion[seriation[left.border],]))
    left.val <- criterion[seriation[left.border],left.max]

    ### break ties
    #if(left.val == right.val) left <- runif(1) > 0.5
    #else (left <- left.val > right.val)

    left <- left.val > right.val

    if(left == TRUE) {
      ### place column left
      chosen <- left.max
      left.border <- left.border - 1
      seriation[left.border] <- chosen
      placed <- c(placed, chosen)
      criterion[,chosen] <- NA


    }else{
      ### place column right
      chosen <- right.max
      right.border <- right.border + 1
      seriation[right.border] <- chosen
      placed <- c(placed, chosen)
      criterion[,chosen] <- NA


    }

    #cat(paste(chosen, if(left == TRUE) paste("left", left.val) else paste("right", right.val), "\n"))
  }

  ### remove empty placeholders
  return(seriation[seriation > 0])
}


