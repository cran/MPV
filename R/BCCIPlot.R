BCCIPlot <-
function(data, k1=1, k2=2, h, h2, output=FALSE, g, layout=FALSE, incl.biasplot=FALSE, plotdata = TRUE) {
    xy <- data
    bias.out <- LPBias(xy,k1=k1,k2=k2,h=h,h2=h2)
    x <- bias.out$x
    y <- bias.out$y
    BCCI.UL <- bias.out$BCCI.UL
    BCCI.LL <- bias.out$BCCI.LL
    biasCI.LL <- bias.out$biasCI.LL
    biasCI.UL <- bias.out$biasCI.UL
    if (layout) {
        oldpar <- par(mfrow=c(2,1), mar=c(4,4,0,1))
        on.exit(par(oldpar))
    }
    y.range <- range(c(BCCI.UL, BCCI.LL))
    plot(x, BCCI.LL, ylim=y.range, type="l", ylab=names(xy)[2], xlab=
names(xy)[1], lty=2)
    lines(x,bias.out$yBC)
    lines(x, BCCI.UL, lty=2)
    if (plotdata) points(data[,1], data[,2])
    if (!missing(g)) lines(x, g(x), col=2, lwd=2)
    if (incl.biasplot) {
        plot(x, biasCI.LL, ylim=y.range-mean(y), type="l", xlab=names(xy)[1], 
ylab="bias")
        lines(x, biasCI.UL)
    }
    Bias <- function(x,y,k,numgrid=401,g, h){
         n <- length(x);     a <- min(x);    b <- max(x); # h <- dpill(x,y)
         g0.l1 <- locpoly(x, g(x), bandwidth=h, degree=k, range.x=c(a-h,b+h), gridsize=numgrid)
         bias <-  approxfun(g0.l1$x, g0.l1$y)(g0.l1$x) - g(g0.l1$x)
         list(x=g0.l1$x,y=bias)
    }
    if (!missing(g)) lines(Bias(xy$x,xy$y,k=k1,numgrid=401,g=g,h=h), col=2, lwd=2)
    abline(0,0)
    if (output) {
        bias.out
    } else {
        invisible()
    }
}
