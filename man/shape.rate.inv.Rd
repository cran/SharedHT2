\name{shape.rate.inv}
\alias{shape.rate.inv}
\title{Link function for EB.Anova/general model}
\description{Computes the model coefficients for the EB.Anova/general
  model given the \code{shape} parameter and \code{rate} matrix parameter}
\usage{
  shape.rate.inv(object)
}
\arguments{
  \item{object}{a list containing two components, \code{shape}, a
    scalar, and \code{rate}, the \code{d} by \code{d} rate matrix
    parameter of the inverse Wishart distribution}
}
 
\value{
  A named numeric vector of length \code{d*(d+1)/2 + 1} containing the
  model coefficients for the EB.Anova/general model.
}

\keyword{utilities}
