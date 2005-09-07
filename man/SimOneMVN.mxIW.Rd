\name{SimOneMVN.mxIW}
\alias{SimOneMVN.mxIW}
\title{Simulate a micro-array experiment dataset}
\description{
  \code{SimOneMVN.mxIW} generates a single simulated micro-array expression 
  experiment under the Multivariate Normal/mixed Inverse Wishart model. This
  can be used to generate a variety of example datasets which violate the 
  assumptions of the Multivariate Normal/Inverse Wishart model.  See the 
  documentation for \code{SimMVN.mxIW} for more details (type ?SimMVN.mxIW).}
}

\usage{
  SimOneMVN.mxIW(shape = NULL, rate = NULL, theta = NULL, f1f2 = c(1/4, 1/2),
                 nreps, Ngenes, effect.size)
}
\arguments{
  \item{shape}{The shape parameter for the Inverse Wishart distribution with equal
    averaged-over-genes group variance.}
  \item{rate}{The rate parameter matrix, of dimension \code{d} by \code{d}
    where \code{d} is the number of experimental groups.}
  \item{theta}{Alternatively to specifying \code{shape} and \code{rate}
    above, the user can directly specify the model parameters from which
    \code{shape} and \code{rate} are computed.  Type ?EBfit for more details.}
  \item{nreps}{Number of replicates per group.}
  \item{Ngenes}{Number of rows (or genes) in the dataset (micro-array experiment)}
  \item{effect.size}{A vector of length \code{Ngenes} giving the effect size. 
    Rows with population mean zero (not differentially expressed) are set to zero
    while rows with non-zero population mean (differentially expressed) are set to 
    some non-zero value. For a feeling of corresponding power in the naive F test
    of all means identically zero see the documetation on \code{find.ncp} by typing
    ?find.ncp.}
  \item{f1f2}{A vector of length 2 containing strictly increasing values between 
    0 and 1.  These are additional parameters specifying the mixing proportion for
    the discrete two component mixture and the proportion of the fractional part
    of \code{shape/(2*d+2)} specifying the amount by which one of the two shape
    parameters should exceed 1.  The other shape parameter is derived so that
    the averaged-over-genes group variances are equal to their values under the
    corresponding Multivariate Normal/Inverse Wishart model.}
}

\value{
  A dataframe having \code{Ngenes} rows and \code{nreps * d} columns where \code{d}
  is implicit in the dimension of Lambda, (see above).  See the documentation for
  \code{SimAffyDat} for more details.}
}

\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{EB.Anova}}, \code{\link{EBfit}}, \code{\link{SimAffyDat}},
         \code{\link{TopGenes}}, \code{\link{SimNorm.IG}}, \code{\link{SimMVN.IW}},
         \code{\link{SimMVN.mxIW}}, \code{\link{SimOneNorm.IG}}, 
         \code{\link{SimOneMVN.IW}} 
}

\examples{
  shape <- 9.107182
  rate <- matrix(c(0.12789434, 0.08468535, 0.08468535, 0.12390469), 2, 2)
  Ngenes <- 12625
  nreps <- 3
  nTP <- 100
  effect.size <- c(rep(4.33, nTP), rep(0, Ngenes-nTP))

  MyDat <- SimOneMVN.mxIW(shape=shape, rate=rate, Ngenes=Ngenes, nreps=nreps, 
                          effect.size=effect.size)

# notice the names given to the columns by default:

  names(MyDat)

# Now try out 'EB.Anova' on your dataset

  fit.MyDat <- EB.Anova(data=MyDat, labels= "log2.grp" \%,\% (1:2), H0="zero.means")

# View the sorted genelist

  TopGenes(fit.MyDat, FDR=0.05, allsig=TRUE)

}

\keyword{datagen}
