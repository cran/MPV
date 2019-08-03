LPBias <-
function(xy,k1,k2,h,h2,numgrid=401,alpha=.95){
  z <- (1-alpha)/2
    x <- xy[,1]
    y <- xy[,2]
    y <- y[order(x)];    x <- sort(x)
    n <- length(x);     a <- min(x);    b <- max(x)
    h1 <- h
    y.l1 <- locpoly(x, y, bandwidth=h1, degree=k1, range.x=c(a-h1,b+h1), 
gridsize=numgrid)  # estimate of interest
    y.l2 <- locpoly(x,y,bandwidth=h2,range.x=c(a-h1,b+h1), 
degree=k2, gridsize=numgrid) # estimate used to estimate bias
    g2hat <- approxfun(y.l2$x,y.l2$y) 
    g1hat <- approxfun(y.l1$x,y.l1$y)
    ghat.l1 <- locpoly(x, g2hat(x), bandwidth=h1, degree=k1, 
range.x=c(a-h1,b+h1), gridsize=numgrid)
    bias.hat <- approxfun(ghat.l1$x, ghat.l1$y)(x) - g2hat(x)
    y <- y[order(x)]
    x <- sort(x)
    y2 <- y[seq(2,length(y),2)]
    y1 <- y[seq(1,length(y)-1,2)]
    y12bar <- (y1+y2)/2
    sigma2 <- mean((y1 - y12bar)^2 + (y2 - y12bar)^2)
    Kx <- function(gridpoint, h, data) diag(dnorm(data-gridpoint, sd=h)) 
    Xx <- function(gridpoint, h, k, data) outer(data-gridpoint, seq(0,k), function(x,y) x^y)
    Cmat <- function(n, k, h, x) {
        C <- matrix(0, nrow=n, ncol=n)
        e1 <- c(1, rep(0,k))
        for (j in 1:n) {
            gridpoint <- x[j]
            X <- Xx(gridpoint, h, k, x)
            K <- Kx(gridpoint, h, x)
            KX <- K%*%X
            XKX <- t(X)%*%KX
            C[,j] <- t(solve(XKX, e1))%*%t(KX)
        }
    C
    }
    C1 <- Cmat(n, k=k1, h=h1, x=x)
    C2 <- Cmat(n, k=k2, h=h2, x=x)
    variance <- apply((C1%*%C2-C2)^2,1,sum)
    bias.var <- variance*sigma2
    bias.sd <- sqrt(bias.var)
    ci.ul <- bias.hat + bias.sd*qnorm(1-z)
    ci.ll <- bias.hat + bias.sd*qnorm(z)
    var.hat <- sigma2*apply(C1, 2, function(x) sum(x^2))    
    sd.hat <- sqrt(var.hat)
    list(biasCI.LL = ci.ll, biasCI.UL =ci.ul, biasEst = bias.hat, bias.sd=bias.sd,
x=x, y=g1hat(x), sd.hat=sd.hat, basicCI.LL = g1hat(x) + sd.hat*qnorm(z), 
basicCI.UL = g1hat(x) + sd.hat*qnorm(1-z))
}
