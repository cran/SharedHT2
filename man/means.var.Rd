\name{means.var}
\alias{means.var}
\title{Means and Variance}
\description{Calculate the vector of group means and empirical
 variance/covariance matrix}
\usage{
 means.var(data, labels, subset, H0 = NULL, Var.Struct = NULL,
                 na.action = na.pass)
}
\arguments{
  \item{data}{data frame containing the logged (base 2 or
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
  \item{subset}{an index vector indicating which rows should be used.
    (NOTE: If given, this argument must be named.)}
  \item{H0}{Character string specifying the null hypothesis. See the
    documentation for \code{EB.Anova}. This argument is tested to
    determine the minimum number of replicates neccesary and is used
    internally in order that rows deleted due to too many missing
    replicates during a call to EB.Anova are also deleted during an
    embedded call to this function. The user can either specify this
    argument or leave it unspecified.}
  \item{Var.Struct}{Character string specifying the variance
    structure. See the documentation for \code{EB.Anova}. This argument
    is tested to determine the minimum number of replicates neccesary
    and is used internally in order that rows deleted due to too many
    missing replicates during a call to \code{EB.Anova} are also deleted
    during an embedded call to this function. The user can either
    specify this argument or leave it unspecified.} 
  \item{na.action}{Specify \code{na.action} = \code{na.pass} if your
    data contains \code{NA}'s and you wish to treat them as missing at
    random. Otherwise, leave it unspecified. See the documentation for
    \code{EB.Anova}.} 
}
 
\value{
  A list containing two components
  \item{mean}{an \code{Ngenes} by \code{d} matrix representing
    the per gene group means}
  \item{var}{a \code{Ngenes} by \code{d*d} matrix representing
    the per gene empirical variance covariance matrices}
}
\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{EB.Anova}}, \code{\link{EBfit}}, \code{\link{SimAffyDat}},
         \code{\link{TopGenes}}, \code{\link{SimNorm.IG}},
         \code{\link{SimMVN.IW}}, \code{\link{SimMVN.mxIW}}, 
         \code{\link{SimOneNorm.IG}}, \code{\link{SimOneMVN.IW}},
         \code{\link{SimOneMVN.mxIW}}
}

\note{Under the Wishart/Inverse Wishart Bayesian model, the expected
  value of the random per gene covariance matrix is equal to
  \code{rate}/(\code{shape} - 2*d - 2). Thus as a consistency check
  you can check the observed mean against the theoretical mean
  as in the following example.
}
\examples{
  data(SimAffyDat)
  fit.SimAffyDat <- EB.Anova(data=SimAffyDat, labels=c("log2.grp" \%,\% (1:2)),
                             H0="zero.means", Var.Struct = "general")
  mv.SimAffyDat <- EBfit(fit.SimAffyDat)
  mv.SimAffyDat$call[[1]] <- as.name("means.var")
  mv.SimAffyDat <- update(mv.SimAffyDat)

  d <- dim(mv.SimAffyDat$mean)[2]

  apply(mv.SimAffyDat$var, 2, FUN=mean)
  c(EBfit(fit.SimAffyDat)$rate/(EBfit(fit.SimAffyDat)$shape - 2*d - 2))
}

\keyword{utilities}
