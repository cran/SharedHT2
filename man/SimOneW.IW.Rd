\name{SimOneW.IW}
\alias{SimOneW.IW}
\title{Simulate a micro-array experiment dataset}
\description{
  \code{SimOneW.IW} generates a single simulated micro-array expression 
  experiment under the Wishart/Inverse Wishart model. This can be used
  to generate a variety of example datasets.}
}
\usage{
SimOneW.IW(nu = NULL, Lambda = NULL, theta = NULL, nreps, Ngenes, 
           effect.size) 
}

\arguments{
  \item{nu}{The shape parameter for the Inverse Wishart distribution}
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

  MyDat <- SimOneW.IW(nu=nu, Lambda=Lambda, Ngenes=Ngenes, nreps=nreps, 
                      effect.size=c(rep(7.5/nreps^0.5, nTP), rep(0, Ngenes-nTP)))

# notice the names given to the columns by default:

  names(MyDat)

# Now try out 'EB.Anova' on your dataset

  fit.MyDat <- EB.Anova(data=MyDat, labels= "log2.grp" \%,\% (1:2), H0="zero.means")

# View the sorted genelist

  TopGenes(fit.MyDat, FDR=0.05, allsig=TRUE)

}

\keyword{datagen}
