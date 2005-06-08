\name{find.ncp}
\alias{find.ncp}
\title{Find non-centrality parameter}
\description{In order to provide a yardstick in designing simulation 
experiments this function calculates the non-centrality parameter for a 
simple F test of the null hypothesis of zero group means}
\usage{
find.ncp(e.I, e.II, nreps, d)
}
\arguments{
  \item{e.I}{type I error i.e. level of significance}
  \item{e.II}{type II error i.e. 1-power}
  \item{nreps}{\# replicates}
  \item{d}{\# groups}
}
 
\value{
  An numeric vector of length 2 containing the calculated non-centrality 
  parameter and the corresponding power of the standard F test of the null
  hypothesis that all group means are zero.
  \item{theta}{The non-centrality parameter}
  \item{power}{Corresponding power}
}
\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{EBfit}}, \code{\link{SimAffyDat}} \code{\link{TopGenes}}, 
         \code{\link{SimW.IW}}, \code{\link{Simnu.mix}}
}

\note{
 Once you derive the non-centrality parameter then you can use it
 in specifying the population means for genes having differential expression.
 For example if investiagation using \code{find.ncp} turns up the value
 7.5 then use $7.5/nreps^0.5$ to specify group means corresponding to 
 experimental conditions under which the gene is differentially expressed.}
}
\examples{
  
  find.ncp(0.005, 0.0316, 3, 2)

}

\keyword{utilities}
\keyword{design}
