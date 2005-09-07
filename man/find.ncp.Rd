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
  An numeric vector of length 2, containing the calculated
  \code{effect-size} and the corresponding power of the standard F test
  of the null hypothesis that all group means are zero at the alternative
  that each group has mean \code{effect-size} times the common within
  group standard deviation.
  \item{theta}{The effect size}
  \item{power}{Corresponding power}
}
\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{SimNorm.IG}}, \code{\link{SimMVN.IW}},
  \code{\link{SimMVN.mxIW}}  
}

\note{
  This is used to set up the true positives genes in the simulation
  routines by deciding which genes are ``true positives" in tests of the null
  that all groups means are zero. Each group is assigned a mean equal to
  $\sigma \theta$ where $\sigma$ is the average common within group
  standard deviation and $\theta$ is the effect size returned by a call
  to this function, for example $\theta=2.0$ as in the example below.
  In calls to the simulation routine you only need to specify the per
  gene effect size as $\sigma$ is calculated internally in all cases.
  For example if you want the first 100 rows to be ``true-positives" in
  simulated micro-array data containing 15000 genes then in the
  simulation routines just specify ``effect.size=c(rep(es, 100), rep(0,
  15000-100)" in the calls to simulation routines listed above. See the
  documentation for routines listed above in ``see also".
}

\examples{
  
  find.ncp(0.005, 0.13179983, 3, 2)

}

\keyword{utilities}
\keyword{design}
