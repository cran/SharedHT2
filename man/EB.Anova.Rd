\name{EB.Anova}
\alias{EB.Anova}
\title{Per Gene Group Comparison with Empirical Bayes Anova Tests}
\description{
  \code{EB.Anova} calculates per gene empirical Bayes Anova
  statistic for testing a variety hypotheses in replicated array
  experiments involving group comparisons (Izmirlian \& Xu, see
  the manuscript in the ./doc directory).  Also appropriate for time
  course data. The Empirical Bayes statistic is calculated using the
  ordinary group means and a shrinkage variance estimate formed as the 
  posterior mean of the variance. If \code{Var.Struc} is set to
  \code{"general"} then this is the posterior mean of the variance/covariance
  matrix under the Wishart/Inverse Wishart Bayesian model. If
  \code{Var.Struct} is set to \code{"simple"} then this is the posterior
  mean of the normalized within group sum of squares under a
  Chi-squared/Inverse Gamma Bayesian model.  In both cases, the
  paremeters of the prior distribution are fit using MLE applied to
  the per gene residual squared error (scalar or matrix corresponding
  to the above).
}
\usage{

  EB.Anova(data, labels, H0 = "equal.means", Var.Struct = "general",
           verbose = T, subset, theta0 = NULL, gradient = F, fit.only = F)
}
\arguments{
  \item{data}{a data frame containing the logged (base 2 or
    base 10) expression values in all arrays from the experiment.
    By default the variables are taken from the environment which
    \code{EB.Anova} is called from.  Variable names should be
    chosen to be internally consistent in some searchable way.
    For example, if you have \code{d}=2 experimental groups (say
    treatment one versus control and treatment two versus control),
    and \code{n}=3 replicates in each group, you might choose names
    like: \code{log2.grp1.n1}, \code{log2.grp2.n1}, \code{log2.grp1.n2},
    \code{log2.grp2.n2}, \code{log2.grp1.n3}, \code{log2.grp2.n3}
    Notice that order that the names occur is irrelevent. In time
    course data the time point is the grouping variable. The rows
    should be named using the gene identifiers.}
  \item{labels}{A character vector containing the group names, these
    being fragments of the variable names in the \code{data} argument
    supplied. In the example above, \code{labels =
    c("log2.grp1", "log2.grp2")}}
  \item{H0}{a character string giving the form of the null hypothesis
    that is to be tested. Specifically
    (i)   if \code{n > max(d, 1+ I(Var.Struct=="general"))}, you can test
    against the \code{H0="zero.means"} null,
    (ii)  if \code{n >= max(d, 1)}, you can test against the
    \code{H0="equal.means"} null or \code{H0="no.trend"} null,
    but of course the latter only makes sense if \code{d > 2}.  Set to
    \code{"equal.means"} by default}
  \item{Var.Struct}{set to either \code{"general"} or \code{"simple"}.
    The default, \code{"general"}, fits the Wishart/Inverse Wishart model
    and computes per gene Empirical Bayes Hotelling T-Squared tests.
    The \code{"simple"} option assumes equal group variances, fits the
    Chi-Squared/Inverse Gamma model and computes per gene Empirical
    Bayes F-tests (or Univariate T-Squared)}
  \item{verbose}{Do you want a trace of the optimization procedure. Set
    to \code{TRUE} by default.}
  \item{subset}{an index vector indicating which rows should be used.
    (NOTE: If given, this argument must be named.)}
  \item{theta0}{optional values for the starting parameters. Must be of
    length \code{d*(d+1)/2 + 1}}
  \item{gradient}{set to \code{TRUE} to optimize using methods requiring
    analytic derivatives. Set to \code{FALSE} by default. The Nelder
    Meade converges from any starting position and fits in less than
    ten seconds on a pentium 4. The likelihood surface is a giant spike
    since there is such an abundance of data on the error structure.}
  \item{fit.only}{set to \code{TRUE} if you only want the result of the
    model fit and not the list of per gene statistics. Set to
    \code{FALSE} by default.}
}

\value{
  An object of class \code{fit.n.data} containing two components:
  \item{data}{A data frame containing the per gene test statistics
    (both the empirical Bayes and the standard versions) together
    with corresponding p-values.}
  \item{EBfit}{An object of class \code{EBfit} containing the results
    of the fitted model. Type \code{?EBfit} for details.}
}

\note{
  Under the model assumptions, the test statistic has an F distribution
  with \code{r} degrees in top and \code{nu + n - 2*r - 2} degrees in bottom.
  where \code{r} = \code{d} for the \code{"zero.means"} test, 
  \code{r} = \code{d-1} for the \code{"equal.means"} test, and
  \code{r} = \code{1} for the \code{"no.trend"} test.
  The test performs quite well, even the asymptotic p-values make
  sense and are 'FDR'-able under a variety of departures from the model.
  It is entirely coded in C.  In an experiment with \code{N} = 12625 genes,
  \code{d}=2 groups with \code{n}=3 replicates per group, the model was
  fit and the list of statistics was computed in less than 10 seconds on a
  pentium 4.
}

\references{
  Izmirlian, G and Xu, J.-L. (2002), \emph{The Shrinkage Variance Hotelling
    T-Squared Test for Genomic Profiling Studies}, NCI technical report.
}
\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{EBfit}}, \code{\link{SimAffyDat}}, \code{\link{TopGenes}}, 
         \code{\link{SimW.IW}}, \code{\link{Simnu.mix}}
}

\examples{

# The included example dataset is a simulated Affymetrix oligonucleotide
# array experiment. Type ?SimAffyDat for details.

  data(SimAffyDat)


# Fit the Wishart/Inverse Wishart empirical Bayes model and derive per gene
# Shared Variance Hotelling T-Squared (ShHT2) statistics.

  fit.SimAffyDat <- EB.Anova(data=SimAffyDat, labels=c("log2.grp" \%,\% (1:2)),
                             H0="zero.means", Var.Struct = "general")

# Top 20 genes (sorted by decreasing ShHT2 statistic) and model summary

  fit.SimAffyDat

# Same screen output & opens html browser with genelist linked to GeneCards database.
# Type ?TopGenes for help

  print(fit.SimAffyDat, browse = TRUE)

# Only the genes selected by the Benjamini-Hochberg procedure at FDR=0.05

  print(fit.SimAffyDat, FDR=0.05, allsig=TRUE)

# Just the top 35 genes

  print(fit.SimAffyDat, n.g = 35)

# In the previous two cases, supplying the argument 'browse'=TRUE produces
# the expected result. If just the genelist without the model summary is desired
# then use calls to 'TopGenes' instead of calls to 'print' in the above
# with exactly the same sytax otherwise.

# Try the update method with Var.Struct="simple":

  fitSV.SimAffyDat <- update(fit.SimAffyDat, Var.Struct = "simple")

# If for some reason, you want the 'EBfit' component then use either

  x <- fit.SimAffyDat$EBfit   # or
  x <- EBfit(fit.SimAffyDat)

# The 'EBfit' print method supplies the simple summary table of coefficients,
# standard errors and Wald statistic p-values mentioned above.

  x

# Notice that the actual structure is more detailed:

  names(x)
  x$log.likelihood
  x$variance
 
# You can perform non-assignment operations directly on the data component
# of the object using 'as.data.frame', which contains the statistics,
# unsorted in the  same order as the original dataset, using 'as.data.frame':

  as.data.frame(fit.SimAffyDat)[1:100, ]

}
\keyword{models}
\keyword{multivariate}
\keyword{univar}
\keyword{robust}
\keyword{htest}
