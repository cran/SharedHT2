\name{SimMVN.IW}
\alias{SimMVN.IW}
\title{Simulation study using the Multivariate Normal/Inverse Wishart Model}
\description{
  \code{SimMVN.IW} generates \code{nsim} simulation replicates of a 
  simulated micro-array expression experiment under the Multivariate Normal/
  Inverse Wishart model for the purposes of benchmarking the empirical Bayes 
  Hotelling T-squared test against three other statistics, the standard
  Hotelling T-squared test, the empirical Bayes univariate F test and the
  standard univariate F test.}
}
\usage{
SimMVN.IW(nsim, shape = NULL, rate = NULL, theta = NULL, nreps, 
          Ngenes, effect.size, FDRlist = 0.05 * (1:5), verbose = F, 
          gradient = F) 
}

\arguments{
  \item{nsim}{Number of simulations desired.}
  \item{shape}{The shape parameter for the Inverse Wishart distribution}
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
  \item{FDRlist}{A list of FDR values to use in the BH step down procedure
    used in summarizing the simulation results.}
  \item{verbose}{Doesn't really belong here. Defaults to FALSE, leave it that
    way or your log file will grow to large}
  \item{gradient}{Also doesn't belong here. Defaults to FALSE, leave it that 
    way.}
}

\value{
  A list containing the following 5 components:
  \item{fdrtbl}{A matrix of dimension \code{length(FDRlist)} by 8, with one row
    corresponding to each of the expected FDR's given in FDRlist and having columns
    \code{ShHT2-TPR}, \code{ShHT2-FPR}, \code{HT2-TPR}, \code{HT2-FPR}, 
    \code{ShUT2-TPR}, \code{ShUT2-FPR},  \code{UT2-TPR}, \code{UT2-FPR}, each giving
    the average empirical true/false positive rate over the \code{nsim} simulations 
    for the corresponding statistic.}
  \item{countstbl}{An \code{Ngenes} by 8 matrix. Instead of thresh-holding the 
    corresponding p-values by the BH stepdown criterion, unique values of the 
    statistic are treated as candidate threshold values giving at each simulation 
    rep empirical true/false positive rates for each of the 4 statistics. These 
    values are averaged over the \code{nsim} simulation reps producing an \code{Ngenes} 
    by 8 matrix.}
  \item{coef}{An \code{nsim} by \code{d*(d+1)/2 + 1} matrix containing the fitted
    model coefficients for the Multivariate Normal/Inverse Wishart model.}
  \item{coefEV}{An \code{nsim} by 2 matrix containing the fitted model coefficients
    for the Normal/Inverse Gamma model.}
  \item{call}{The original call to \code{SimMVN.IW}}
}

\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{EB.Anova}}, \code{\link{EBfit}}, \code{\link{SimAffyDat}},
         \code{\link{TopGenes}}, \code{\link{SimNorm.IG}}, \code{\link{SimMVN.mxIW}}, 
         \code{\link{SimOneNorm.IG}}, \code{\link{SimOneMVN.IW}},
         \code{\link{SimOneMVN.mxIW}}
}

\examples{\dontrun{
  data(SimAffyDat)
  fit.SimAffyDat <- EB.Anova(data=SimAffyDat, labels=c("log2.grp" %,% (1:2)),
                             H0="zero.means", Var.Struct = "general")

  SimResults <- 
    SimMVN.IW(nsim=500, theta=EBfit(fit.SimAffyDat)$coef, Ngenes=12625,
              nreps=3, FDRlist = 0.05*(1:5), effect.size = c(rep(4.33, 100), 
              rep(0, 12625 - 100)))
          
# Or create a batch file like this
# contents of mysim.R:
  library(SharedHT2)
  nsim <- 500
  nreps <- 3
  Ngenes <- 12625
  nTP <- 100
  effect.size <- c(rep(4.33, nTP), rep(0, Ngenes - nTP)
  theta <- EBfit(fit.SimAffyDat)$coef

  SimResults <- SimMVN.IW(nsim=nsim, theta = theta, nreps = nreps, 
                          Ngenes = Ngenes, effect.size = effect.size)

# At the command prompt

  R CMD BATCH mysim.R mysim.Rout

# nsim=500, Ngenes=12625, nreps=3, with d=2 groups (implicit in the dimension
# of theta) will take just under 3 hours on a pentium 4.
}}

\keyword{datagen}
