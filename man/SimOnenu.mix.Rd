\name{SimOnenu.mix}
\alias{SimOnenu.mix}
\title{Simulate a micro-array experiment dataset}
\description{
  \code{SimOnenu.mix} generates a single simulated micro-array expression 
  experiment under the Wishart/mixed Inverse Wishart model. This can be used
  to generate a variety of example datasets which violate the assumptions
  of the Wishart/Inverse Wishart model.  See the documentation for 
  \code{Simnu.mix} for more details (type ?Simnu.mix).}
}

\usage{
  SimOnenu.mix(nu = NULL, Lambda = NULL, theta = NULL, f1f2 = c(1/4, 1/2),
               nreps, Ngenes, effect.size)
}
\arguments{
  \item{nu}{The shape parameter for the Inverse Wishart distribution with equal
    averaged-over-genes group variance.}
  \item{Lambda}{The rate parameter matrix, of dimension \code{d} by \code{d}
    where \code{d} is the number of experimental groups.}
  \item{theta}{Alternatively to specifying \code{nu} and \code{Lambda}
    above, the user can directly specify the model parameters from which
    \code{nu} and \code{Lambda} are computed.  Type ?EBfit for more details.}
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
    of \code{nu/(2*d+2)} specifying the amount by which one of the two shape
    parameters should exceed 1.  The other shape parameter is derived so that
    the averaged-over-genes group variances are equal to their values under the
    corresponding Wishart/Inverse Wishart model.}
}

\value{
  A dataframe having \code{Ngenes} rows and \code{nreps * d} columns where \code{d}
  is implicit in the dimension of Lambda, (see above).  See the documentation for
  \code{SimAffyDat} for more details.}
}

\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{EB.Anova}}, \code{\link{EBfit}}, \code{\link{SimAffyDat}},
         \code{\link{TopGenes}}, \code{\link{Simnu.mix}}
}

\examples{
  nu <- 9.107182
  Lambda <- matrix(c(0.12789434, 0.08468535, 0.08468535, 0.12390469), 2, 2)
  Ngenes <- 12625
  nreps <- 3
  nTP <- 100

  MyDat <- SimOnenu.mix(nu=nu, Lambda=Lambda, Ngenes=Ngenes, nreps=nreps, 
                        effect.size=c(rep(7.5/nreps^{0.5}, nTP), rep(0, Ngenes-nTP)))

# notice the names given to the columns by default:

  names(MyDat)

# Now try out 'EB.Anova' on your dataset

  fit.MyDat <- EB.Anova(data=MyDat, labels= "log2.grp" \%,\% (1:2), H0="zero.means")

# View the sorted genelist

  TopGenes(fit.MyDat, FDR=0.05, allsig=TRUE)

}

\keyword{datagen}
