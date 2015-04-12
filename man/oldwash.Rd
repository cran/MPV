\name{oldwash}
\alias{oldwash}
\docType{data}
\title{oldwash}

\description{
The \code{oldwash} dataframe has 49 rows and 8 columns. 
The data are from the start up of a wash still considering the amount of time it takes to heat up to a specified temperature and possible influencing factors. 
}
\usage{data("oldwash")}
\format{
  A data frame with 49 observations on the following 8 variables.
  \describe{
    \item{\code{Date}}{character, the date of the run}
    \item{\code{startT}}{degrees Celsius, numeric, initial temperature}
    \item{\code{endT}}{degrees Celsius, numeric, final temperature}
    \item{\code{time}}{in minutes, numeric, amount of time to reach final temperature}
    \item{\code{Vol}}{in litres, numeric, amount of liqiud in the tank (max 2000L)}
    \item{\code{alc}}{numeric, the percentage of alcohol present in the liquid}
    \item{\code{who}}{character, relates to the person who ran the still}
    \item{\code{batch}}{factor with levels 1 = first time through, 2 = second time through}
  }
}
\details{
The purpose of the wash still is to increase the percentage of alcohol and strip out unwanted particulate. 
It can take a long time to heat up and this can lead to problems in meeting production time limits. 
}
\source{
Charisse Woods, Endless Summer Distillery (2014)
}
\examples{

oldwash.lm<-lm(log(time)~startT+endT+Vol+alc+who+batch,data=oldwash)
summary(oldwash.lm)
par(mfrow=c(2,2))
plot(oldwash.lm)

data2<-subset(oldwash,batch==2)
hist(data2$time)
data1<-subset(oldwash,batch==1)
hist(data1$time)

oldwash.lmc<-lm(time~startT+endT+Vol+alc+who+batch,data=data1)
summary(oldwash.lmc)
plot(oldwash.lmc)

oldwash.lmd<-lm(time~startT+endT+Vol+alc+who+batch,data=data2)
summary(oldwash.lmd)
plot(oldwash.lmd)

}
\keyword{datasets}
