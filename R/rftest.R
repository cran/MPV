rftest <- function(u, m=5) {
    n <- length(u)-1
    A <- diag(rep(1,n))
    A <- rbind(rep(0, n), A)
    A <- cbind(A, rep(0, n+1))
    xy <- matrix(u, nrow=n+1)
    for (j in 1:m) {
        xy <- cbind(xy, A%*%xy[,j])
    }
    xy <- data.frame(xy)
    names(xy) <- c("y", paste("x", 1:m, sep=""))
    xy <- xy[-(1:m),]
    xytrain <- xy[1:(n/2),]
    xytest <- xy[-(1:(n/2)),]
    xy.rf <- randomForest(y ~ ., data = xytrain)
    oldpar <- par(mfrow=c(1,2))
    on.exit(par(oldpar))
    plot(predict(xy.rf), xytrain$y, cex=.3, 
    xlab="predicted values", ylab="observed values",
        main = "training data")
    abline(lm(xytrain$y ~predict(xy.rf)))
    plot(predict(xy.rf, newdata=xytest), xytest$y, cex=.3,
    xlab="predicted values", ylab="observed values", 
        main="test data")
    abline(lm(xytest$y ~ predict(xy.rf, newdata=xytest)))
}
