
# coding.r
#
# dummy coding for data mining applications
#
# fixme: no reverse methods implemented
#
# ceeboo 2005

as.dummy <- function(x, ...) 
    UseMethod("as.dummy")

as.dummy.logical <- function(x, ...)
    x

as.dummy.integer <- function(x, ...) { 
    x <- as.dummy(as.factor(x))
    x
}

as.dummy.factor <- function(x, ...) {
    x <- .Call("asDummy", x)
    x
}

as.dummy.matrix <- function(x, sep=" ", ...) {
    if (is.null(colnames(x)))
       colnames(x) <- paste("V", 1:dim(x)[2], sep="")
    obj <- NULL
    colnames <- NULL
    levels <- NULL
    for (i in 1:dim(x)[2]) {
        z <- as.dummy(x[,i])
        obj <- cbind(obj, z)
        levels <- c(levels, list(levels(z)))
        colnames <- c(colnames, paste(colnames(x)[i], levels(z), sep=sep))
    }
    rownames(obj) <- rownames(x)
    colnames(obj) <- colnames
    names(levels) <- colnames(x)
    attr(obj, "levels") <- levels
    obj
}

as.dummy.list <- function(x, ...)
    lapply(x, function(z) as.dummy(z))

as.dummy.data.frame <- function(x, sep=" ", ...) {
    if (is.null(names(x)))
       names(x) <- paste("V", 1:length(x), sep="")
    obj <- NULL
    levels <- NULL
    colnames <- NULL
    for (name in names(x)) {
        z <- as.dummy(x[[name]])
        obj <- cbind(obj, z)
        levels <- c(levels, list(levels(z)))
        colnames <- c(colnames, paste(name, levels(z), sep=sep))
    }
    if (is.data.frame(x))
       rownames(obj) <- rownames(x)
    colnames(obj) <- colnames
    names(levels) <- colnames(x)
    attr(obj, "levels") <- levels
    obj
}

###
