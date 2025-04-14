\name{postunstack}
\alias{postunstack}
\title{
Unstack Vectors into a Data Frame 
}
\description{
Padding an unstacked data frame with missing values to ensure equal length vectors in resulting list.
This list is then coerced into a data frame for ease of producing tables.  
}
\usage{
postunstack(x, form, ...)
}
\arguments{
  \item{x}{A list or data frame to be stacked or unstacked.}
  \item{form}{a two-sided formula whose left side evaluates to the vector
          to be unstacked and whose right side evaluates to the
          indicator of the groups to create.  Defaults to 'formula(x)'
          in the data frame method for 'unstack'.}
  \item{...}{further arguments passed to or from other methods.}
}
\value{a data frame of columns according to the formula
     'form'.  If the columns do not all have the same length, the resulting
     list is coerced to a data frame by padding with missing values.
}
\author{
W. John Braun
}
\seealso{\code{\link{unstack}}}
\keyword{manip}
