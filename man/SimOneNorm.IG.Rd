\name{SimOneNorm.IG}
\alias{SimOneNorm.IG}
\title{Simulate a micro-array experiment dataset}
\description{
  \code{SimOneNorm.IG} generates a single simulated micro-array expression 
  experiment under the Normal/Inverse Gamma model. This can be used
  to generate a variety of example datasets.}
}
\usage{
SimOneNorm.IG(shape = NULL, rate = NULL, theta = NULL, ngroups, nreps, Ngenes, 
              effect.size) 
}

\arguments{
  \item{shape}{The shape parameter for the Inverse Gamma distribution}
  \item{rate}{The rate parameter for the Inverse Gamma distribution}
  \item{theta}{Alternatively to specifying \code{shape} and \code{rate}
    above, the user can directly specify the model parameters, i.e. the
    logged shape and logged rate.}
  \item{ngroups}{The number of experimental groups.}
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
         \code{\link{TopGenes}}, \code{\link{SimNorm.IG}}, \code{\link{SimMVN.IW}},
         \code{\link{SimMVN.mxIW}}, \code{\link{SimOneMVN.IW}}, 
         \code{\link{SimOneMVN.mxIW}}
}

\examples{
\dontrun{
  shape <- 1.93589032
  rate <- 0.04020591
  Ngenes <- 12625
  ngroups <- 2
  nreps <- 3
  nTP <- 100
  effect.size <- c(rep(4.33, nTP), rep(0, Ngenes-nTP))

  MyDat <- SimOneNorm.IG(shape=shape, rate=rate, ngroups=ngroups, Ngenes=Ngenes, nreps=nreps, 
                         effect.size=effect.size)

# notice the names given to the columns by default:

  names(MyDat)

# Now try out 'EB.Anova' on your dataset

  fit.MyDat <- EB.Anova(data=MyDat, labels= "log2.grp" \%,\% (1:2), Var.Struct="simple",
                        H0="zero.means")

# View the sorted genelist

  TopGenes(fit.MyDat, FDR=0.05, allsig=TRUE)

}
}
\keyword{datagen}
