\name{EBfit}
\alias{EBfit}
\title{Extractor function}
\description{Extract the fit component from an object of class \code{fit.n.data}}
\usage{
  EBfit(object)
}
\arguments{
  \item{object}{An object of class \code{fit.n.data} returned by EB.Anova}
}
 
\value{
  An object of class \code{EBfit} containing the following components:  
  \item{coefficients}{a numerical vector of length \code{d*(d+1)/2 + 1}
    containing the model parameters.  The first is the logged shape
    parameter of the inverse Wishart, and the remaining \code{d*(d+1)/2} are
    related to the its rate matrix parameter, which is symmetric and
    \code{d} by \code{d}.  The first \code{d} of these are the logged
    diagonal elements and the remaining \code{d*(d-1)/2} are the
    off-diagonal elements.}
  \item{variance}{the estimated \code{d*(d+1)/2 + 1} by \code{d*(d+1)/2 + 1}
    covariance matrix of the model coefficients}
  \item{gradient}{the value of the gradient vector at the MLE--hopefully
    very close to zero.}
  \item{nu}{the shape parameter of the inverse Wishart =
    \code{exp(coefficients[1])}}
  \item{Lambda}{the rate parameter of the inverse Wishart, a positive
    definite, symmeteric \code{d} by \code{d} matrix.}
  \item{likelihood}{the log likelihood...actually positive (like I said,
    big spike)}
  \item{fail, fncount, grcount, mask}{details from the optimization
    routine}
  \item{usegr}{whether the gradient method was used in the optimization}
  \item{call}{the original call to \code{EB.Anova}}  
}
\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{EBfit}}, \code{\link{SimAffyDat}} \code{\link{TopGenes}}, 
         \code{\link{SimW.IW}}, \code{\link{Simnu.mix}}
}
  
  \examples{
#  See examples under ?EB.Anova  
}
\keyword{utilities}
