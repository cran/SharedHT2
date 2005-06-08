\name{shape.rate}
\alias{shape.rate}
\title{Inverse link for EB.Anova/general model}
\description{Computes the \code{shape} parameter and \code{rate} matrix
  parameter given the coefficients from an EB.Anova/general, i.e.
  Wishart/Inverse Wishart model}
\usage{
  shape.rate(object)
}
\arguments{
  This function expects an object belonging to any one of the folloiwng
  classes: \code{fit.n.data}, \code{EBfit}, or \code{numeric}. If
  \code{numeric} then its length should be \code{d*(d+1)/2 + 1} for
  some \code{d}
  \item{object}{see above}
}
 
\value{
  An list containing 2 components:
  \item{shape}{a scalar equal to the shape parameter of the inverse
    Wishart distribution. The larger it is the more heteroscedastic 
    the extra-gene variability}
  \item{rate}{the \code{d} by \code{d} rate matrix parameter of the
    inverse Wishart distribution}
}

\keyword{utilities}
