"lof.lm" <-
function (lm.obj, approx = FALSE, call.plot=TRUE) 
{
   Xmatrix <- model.matrix(lm.obj)
   p <- dim(Xmatrix)[2]
   if (p==2) {x <- Xmatrix[,-1]}
       else {x <- Xmatrix[,1]}
   y <- model.response(lm.obj$model, "numeric")
   if (call.plot) {plot(y~x)
                   abline(lm.obj)}
   if (approx) {
      cat("The following results are only approximate!!! \n")
      x.sort <- sort(x)
      y <- y[order(x)] 
      n <- length(x)
      new.indices <- rep(seq(1,(n+1)/2),rep(2,(n+1)/2))[1:n]
      x <- rep(sapply(split(x.sort,new.indices),mean),rep(2,(n+1)/2))[1:n]
      lm.obj <- lm(y~x)
      if (call.plot) {points(x,y,pch=16,col=2)
                      abline(lm.obj,col=2)
                      title(sub="red denotes approximation")}
    }  
   y.x <- split(y, x)
   if (call.plot) {group.means <- sapply(y.x, mean)
                   points(sort(unique(x)),group.means, pch=16, col=4)}
   dfPE <- sum(sapply(y.x, length))-length(y.x)
   if (dfPE == 0) {print("There are no replicate observations.")
                   print("Exact Lack of Fit Test is Not Applicable.")
                   print("For an approximation, try approx = TRUE")
}
   else {
   SS <- function (x) 
         {
          xbar <- mean(x)
          sum((x-xbar)^2)
         }
   
   SSRes <- summary(lm.obj)$sigma^2*(length(y)-p)
   SSPE <- sum(sapply(y.x, SS))
   SSLOF <- SSRes - SSPE
   dfLOF <- length(y) - dfPE - p
   df <- c(dfLOF, dfPE)
   ss <- c(SSLOF, SSPE)
   ms <- ss/df
   f  <- c(ms[1]/ms[2],NA)
   pv <- c(1 - pf(f,df[1],df[2]))
   pred.ratio <- c(diff(range(predict(lm.obj)))/sqrt(p*ms[2]/length(x)),NA)
   table <- data.frame(df, ss, ms, f, pv, pred.ratio)
   dimnames(table) <- list(c("Lack of Fit", "Pure Error"), c("Df", 
        "Sum Sq", "Mean Sq", "F value", "Pr(>F)", "prediction ratio"))
   structure(table, heading = c("Test of Lack of Fit for Simple Linear Regression\n", 
   paste("Response:", deparse(formula(lm.obj)[[2]]))), class = c("anova","data.frame"))
   } 
}
