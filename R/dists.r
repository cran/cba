
# wrapper for various auto- and cross-distance functions.
# mostly compatible with dist.
#
# note that some of the checks may be unecessary, but we 
# cannot rely on R implementation details. what a mess!
#
# fixme: storage mode coerces which may be inappropriate.
#        we should raise an error if this would result in
#        a narrowing conversions, so that the user can
#        decide what she wants.
#
# (C) ceeboo 2005

dists <- function(x, y=NULL, method="minkowski", p=2) {
    METHODS <- c("pdist","pdist","mdist","cdist","bdist","ebdist","fbdist",
                  "adist")
    names(METHODS) <- c("minkowski","manhatten","maximum","canberra","binary",
                        "ebinary","fbinary","angular")
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
# NAs for some entries because we do not have diagonal 
# entries. fixed NA subscripts.
#
# fixme: coercing to double is conceptually not quite 
#        correct but not trivial to implement otherwise.
#
# ceeboo 2005, 2006

subset.dist <- function(x, subset, ...) {
    if (missing(subset))
       return(x)
    # subscript hack
    labels <- array(1:attr(x, "Size"))
    rownames(labels) <- attr(x, "Labels")
    if (!is.double(x))
       storage.mode(x) <- "double"
    # sets Size and Labels
    obj <- .Call("subset_dist", x, subset, labels)
    obj <- structure(obj, class="dist", 
                          Diag=attr(x, "Diag"), Upper=attr(x, "Upper"), 
                          method=attr(x, "method"))
    obj
}

"[[.dist" <- subset.dist

rowSums.dist <- colSums.dist <- function(x, na.rm = FALSE) {
    if (!inherits(x, "dist"))
        stop("'x' not of class 'dist'")
    storage.mode(x) <- "real"
    if (!is.logical(na.rm))
        stop("'na.rm' not logical")
    obj <- .Call("rowSums_dist", x, na.rm)
    names(obj) <- attr(x, "Labels")
    obj
}

dapply <- function(x, y = NULL, FUN, ...) {
    if (!is.matrix(x))
        stop(gettext("'x' not a matrix"))
    storage.mode(x) <- "real"
    if (!is.null(y)) {
        if (!is.matrix(y))
            stop(gettext("'y' not a matrix"))
        if (dim(x)[2] != dim(y)[2])
            stop(gettext("'x' and 'y' do not conform"))
        storage.mode(y) <- "real"
    }
    if (!is.function(FUN))
        stop("'FUN' not a function")
    obj <- .External("apply_dist", x, y, FUN, ...)
    if (!is.null(y)) {
        rownames(obj) <- rownames(x)
        colnames(obj) <- rownames(y)
    } else 
        obj <- structure(obj, Size=dim(x)[1], class="dist",
                              Diag=FALSE, Upper=FALSE,
                              Labels=rownames(x))
    obj
}

##

cluster.dist <- function(x, beta) {
    if (!inherits(x, "dist"))
        stop("'x' not of class dist")
    storage.mode(x) <- storage.mode(beta) <- "real"
    obj <- .Call("cluster_dist", x, beta)
    names(obj) <- attr(x,"Labels")
    obj
}

###

