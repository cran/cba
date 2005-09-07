
# wrapper for various auto- and cross-distance functions.
# mostly compatible with dist.
#
# (C) ceeboo 2005

dists <- function(x, y=NULL, method="minkowski", p=2) {
    METHODS <- c("pdist","pdist","mdist","cdist","bdist","ebdist","adist")
    names(METHODS) <- c("minkowski","manhatten","maximum","canberra","binary",
                        "ebinary","angular")
    if (length(method) != 1 || !method %in% names(METHODS))
       stop("method not implemented")
    if (!is.matrix(x))
       stop(paste(sQuote("x"),"not a matrix"))
    if (method == "binary") {
       if (!is.logical(x))
          stop(paste(sQuote("x"),"not logical"))
    }
    else
       if (!is.real(x))
          storage.mode(x) <- "real"
    if (!is.null(y)) {
       if (!is.matrix(y))
          stop(paste(sQuote("y"),"not a matrix"))
       if (dim(x)[2] != dim(y)[2])
          stop("column dimensions do not conform")
       if (method == "binary") {
          if (!is.logical(y))
             stop(paste(sQuote("y"),"not logical"))
       }
       else
          if (!is.real(y))
             storage.mode(y) <- "real"
    }
    if (method == "minkowski") {
       if (p < 0 || p == Inf)
          stop(paste(sQuote("p"),"illegal value"))
       storage.mode(p) <- "real"

       obj <- .Call(METHODS[method], x, y, p)
       attr(obj,"p") <- p
    }
    else if (method == "manhatten") {
       p <- 0
       storage.mode(p) <- "real"

       obj <- .Call(METHODS[method], x, y, p)
       attr(obj,"p") <- p
    }
    else 
       obj <- .Call(METHODS[method], x, y)
    #
    if (!is.null(y)) {
       rownames(obj) <- rownames(x)
       colnames(obj) <- rownames(y)		  
    }
    else
       obj <- structure(obj, Size=dim(x)[1], class="dist", 
                        Diag=FALSE, Upper=FALSE, 
                        Labels=rownames(x), method=method)
    obj
}

# subsetting for objects of class dist. 
# 
# note that non-unique indexing is allowed but will return
# NAs for some entries because we usually do not have the
# diagonal entries.
#
# ceeboo 2005

subset.dist <- function(x, subset, ...) {
    if (missing(subset))
       return(x)
    if (attr(x, "Diag"))		    # rare, so forget about it
       stop("diagonal not implemented")
    # subscript hack
    labels <- array(1:attr(x, "Size"))
    rownames(labels) <- attr(x, "Labels")
    if (!is.double(x))
       storage.mode(x) <- "double"
    # sets the Labels attribute
    obj <- .Call("subset_dist", x, subset, labels)
    obj <-  structure(obj, Size=length(subset), class="dist", 
                      Diag=FALSE, Upper=attr(x, "Upper"), 
                      method=attr(x, "method"))
    obj
}

"[[.dist" <- subset.dist

###

