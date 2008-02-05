\name{SimAffyDat}
\alias{SimAffyDat}
\title{Simulated Affymetrix(TM) Oligonucleotide Array Experiment Data}
\description{
  The package \code{SharedHT2} contains one dataset \code{SimAffyDat} 
  for testing purposes.  The dataframe has 12625 rows and 6 columns.
  The rows are named after the hgu95 affy gene identifiers in order
  that the links to the GeneCard database work. The six column variables
  represent hypothetical base two logged expression ratios for 
  three replicates of treatment 1 versus control and treatment 2 
  versus control.  
}
\usage{
 data(SimAffyDat)
}
\note{
  All but the first 100 rows have population mean zero.
  The first 100 rows have group specific population means equal to 
  $4.33$ times the square root of the averaged sample variances across 
  genes. The sampling distribution used in the simulation is the
  Wishart/Inverse Wishart model. See the manuscript "sharedHT2.pdf" in
  the ./doc directory. 
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
