\name{SimAffyDat.ann}
\alias{SimAffyDat.ann}
\title{Simulated Affymetrix(TM) Oligonucleotide Array Experiment Data,
  with annotation name attribute.}
\description{
  The package \code{SharedHT2} contains one dataset \code{SimAffyDat} 
  for testing purposes.  This is a copy which illustrates the enhanced
  functionality of the browsing functions for users wishing to install
  the bioconductor packages \pkg{affy} and \pkg{hgu95av2}.  It is an
  identiccal copy of \code{SimAffyDat}, with the only difference being
  an attribute, named "annotation", containing the name of the
  appropriate bioconductor annotation package (in this case
  \pkg{hgu95av2}).  See the help for \code{SimAffyDat}.
}
\usage{
 data(SimAffyDat.ann)
}
\note{
  See the related example under the help in \code{EB.Anova} and
  \code{TopGenes} in order to learn how to use the annotation attribute
  feature.  
}

\references{
  Izmirlian, G and Xu, J.-L. (2002), \emph{The Shrinkage Variance Hotelling
    T-Squared Test for Genomic Profiling Studies}, NCI technical report.
}
\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{EB.Anova}}, \code{\link{EBfit}}, \code{\link{SimAffyDat}},
         \code{\link{TopGenes}}, \code{\link{SimNorm.IG}},
         \code{\link{SimMVN.IW}}, \code{\link{SimMVN.mxIW}}, 
         \code{\link{SimOneNorm.IG}}, \code{\link{SimOneMVN.IW}},
         \code{\link{SimOneMVN.mxIW}}

}

\keyword{datasets}
