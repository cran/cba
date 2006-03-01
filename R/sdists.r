
# implements a wrapper to distance (similarity) computation on 
# collections of sequences. auto and cross distances can be
# computed (compare with dists)
#
# note that 1) we can supply lists of vectors or vectors of 
#              character (strings)
#           2) operation weights are in the order of
#              insertion/deletion, equality, and replacing
#           3) the first row/column of the matrix of alphabet
#              weights are used for replacement with the empty
#              symbol (space)
#           4) include NA if exclude = NULL
#           5) the C function returns NA if NAs are encounterd
#
# ceeboo 2006

sdists <- function(x,y=NULL, method="ow", weight=c(1,0,2),
                   exclude=c(NA,NaN,Inf,-Inf)) {
    METHODS <- c("ow","aw","awl")
    code <- pmatch(method, METHODS)
    if (is.na(code))
       stop("invalid method")
    if (code == -1)
       stop("ambiguous method")
    if (is.character(x))
       x <- strsplit(x,"")
    if (!is.list(x))
       stop("'x' not a list")
    if (code == 2) {
       if (!is.matrix(weight))
          stop("'weight' not a matrix")
       if (dim(weight)[1] != dim(weight)[2])
          stop("'weight' not square")
       if (is.null(colnames(weight)))
          stop("'weight' no colnames")
       l <- colnames(weight)[-1]
    }
    else {
       if (length(weight) < 3)
          stop("'weight' invalid")
       l <- unique(unlist(c(x,y)))          # determine symbol set
    }
    x <- lapply(x,function(x,l) 
             factor(x,levels=l,exclude=exclude), l)
    if (!is.null(y)) {
       if (is.character(y))
          y <- strsplit(y,"")
       if (!is.list(y))
          stop("'y' not a list")
       y <- lapply(y,function(x,l) 
                factor(x,levels=l,exclude=exclude), l)
    }
    storage.mode(weight) <- "real"
    obj <- .Call("sdists",x,y,as.integer(code),weight)
    if (is.null(y))
       obj <- structure(obj, Size=length(x), class="dist",
                             Diag=FALSE, Upper=FALSE,
                             Labels=names(x), method=method)
    else {
       rownames(obj) <- names(x)
       colnames(obj) <- names(y)
    }
    obj
}

###
