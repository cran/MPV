\name{p2.14}
\alias{p2.14}
\non_function{}
\title{Data Set for Problem 2-14}
\usage{data(p2.14)}
\description{
The \code{p2.14} data frame has 8 observations on the molar
ratio of sebacic acid and the intrinsic viscosity of copolyesters.
One is interested in predicting viscosity from the sebacic acid ratio.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{ratio}{molar ratio}
    \item{visc}{viscosity}
  }
}
\source{
 Montgomery, D.C., Peck, E.A., and Vining, C.G. (2001)
Introduction to Linear Regression Analysis. 3rd Edition, John Wiley and Sons.
}
\reference{
Hsuie, Ma, and Tsai (1995) Separation and Characterizations of 
Thermotropic Copolyesters of p-Hydroxybenzoic Acid, Sebacic Acid and
Hydroquinone. Journal of Applied Polymer Science, 56, 471-476.
}
\examples{
data(p2.14)
attach(p2.14)
plot(p2.14, pch=16, ylim=c(0,1))
visc.lm <- lm(visc ~ ratio)
summary(visc.lm)
visc.conf <- predict(visc.lm, interval="confidence")
lines(ratio, visc.conf[,2], col="red")
lines(ratio, visc.conf[,3], col="red")
visc.pred <- predict(visc.lm, interval="prediction")
lines(ratio, visc.pred[,2], col="blue")
lines(ratio, visc.pred[,3], col="blue")
detach(p2.14)
}
\keyword{datasets}
