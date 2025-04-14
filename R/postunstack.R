postunstack <- function(x, form, ...) {
    stkout <- unstack(x, form, ...)
    if (!is.data.frame(stkout)) {
        sizes <- sapply(stkout, length)
        maxsize <- max(sizes)
        for (k in 1:length(sizes)) {
            stkout[[k]] <- c(stkout[[k]], rep(NA, maxsize - sizes[k]))
        } 
        stkout <- as.data.frame(stkout)
    }
    return(stkout)
}
