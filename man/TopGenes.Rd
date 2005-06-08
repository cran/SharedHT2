\name{TopGenes}
\alias{TopGenes}
\title{Gene lists linked to Gene Cards}
\description{
  \code{TopGenes} creates a genelist sorted on values of a chosen statistic.
}
\usage{
  TopGenes(obj, by = "EB", ref = 1, FDR = 0.05, allsig = FALSE, n.g = 20, 
           browse = FALSE, search.url = genecards, path = "", file = "") 
}
\arguments{
  \item{obj}{An object of class \code{fit.n.data} returned by \code{EB.Anova}}
  \item{by}{Specify "EB" or "naive".  Do you want to use the empirical bayes 
    or the naive variant of the statistic. Defaults to "EB".}
  \item{ref}{If you would like to see group to group fold changes included
    in the table then specify \code{ref=k} where \code{k} is the index of the
    group considered referent.}
  \item{FDR}{The false discovery fate (FDR) to use in the Benjamini-Hochberg
    (BH) stepdown procedure.}
  \item{allsig}{Set to \code{TRUE} if you only want to see the genes 
    meeting the BH criterion. Defaults to FALSE}
  \item{n.g}{If \code{allsig} is \code{FALSE}, you must specify the number
    of top genes that you want to see.}
  \item{browse}{set to \code{TRUE} if you the results displayed in the HTML
    browser with gene identifiers linked to the GeneCards database at the 
    Weizmann institute. Defaults to FALSE}
  \item{search.url}{should contain a url href to search an online database
    when a gene identifier is appended onto the end. Defaults to 'genecards',
    included in this package, for searching the GeneCards database at the
    Weizmann institute.}
  \item{file}{If \code{browse} is set to \code{TRUE} and if you want to 
    save the html file then supply a filename, \code{file=foobgwb.html}}
  \item{path}{If specifying a \code{file} name above, optionally you
    may specify a \code{path} name. In the unix implementation, specifying
    a file name without a path writes to the current working directory.
    This feature is not yet supported in Windows. Use \code{path} argument
    to explicitly specify a directory}
}

\value{
  A sorted genelist in the form of a \code{n.g} by \code{5+d} matrix, where \code{d} is
  the number of groups.  The columns contain the following values: 
  following columns:
  \item{RowNum}{Row number from the original unsorted data frame. Useful when you
    happen to know that the first 100 genes are true positives and the rest are not
    (as is the case with the supplied dataset, \code{SimAffyDat}).}
  \item{GeneId}{Taken from the row component of \code{dimnames} in the original
    dataframe.  As such, it is very usefull to name your rows using the affy gene
    identifiers as these are searchable in the GeneCards database.}
  \item{\code{NAME1}}{Group mean corresponding to the group named \code{NAME1}}
  \item{...}{}
  \item{\code{NAMEd}}{Group mean corresponding to the group named \code{NAMEd}}
  \item{\code{TYPE}.stat}{The statistic used in the sort. This is determined
    by user choices at two levels. First, when computations are performed inside
    the call to \code{EB.Anova}, the specification of \code{Var.Struct}, which
    defaults to "general" with alternate value "simple" determins whether the 
    multivariate test, i.e. the Hotelling T-squared, (HT2) or the univariate F
    test (UT2) are computed. In both cases, both the empirical Bayes and standard
    variants are computed. Next after the computation is completed, when the user
    requests the sorted genelist, the specification of the arguement \code{by}, either
    "EB" (default) or "naive" determines which of the two computed statistics is
    used to perform the sort.  Thus \code{TYPE} assumes one of four values
    \code{ShHT2}, \code{HT2}, \code{ShUT2}, or \code{UT2}.  The first two versus
    last two split is on the value of \code{Var.Struct} used in the \code{EB.Anova} 
    computation, while the first/second or third/fourth split is made on the \code{by}
    argument passed to \code{TopGenes}.}
  \item{\code{TYPE}.p-val}{The corresponding p-value under the corresponding model.
    See the manuscript sharedHT2.pdf in the ./doc directory.}
  \item{FDR.stepdown=\code{FDR}}{The BH criterion values computed as 
    \code{FDR * rank/Ngenes}}
}



\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{EB.Anova}}, \code{\link{EBfit}}, \code{\link{SimAffyDat}}, 
         \code{\link{SimW.IW}}, \code{\link{Simnu.mix}}
}

\note{The print method for class \code{fit.n.data} contains a call to 
  \code{TopGenes} so that calls to \code{TopGenes} and to \code{print}
  work almost the same, except the latter produces a model fit summary
  as well.}

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

  TopGenes(fit.SimAffyDat, browse = TRUE)

# Only the genes selected by the Benjamini-Hochberg procedure at FDR=0.05

  TopGenes(fit.SimAffyDat, FDR=0.05, allsig=TRUE)

# Just the top 35 genes

  TopGenes(fit.SimAffyDat, n.g = 35)

# Try the Var.Struct="simple" option:

  fitSV.SimAffyDat <- update(EBfit(fit.SimAffyDat), Var.Struct = "simple")

# Now try TopGenes using the univariate statistic:

  TopGenes(fitSV.SimAffyDat, FDR=0.05, allsig=TRUE)

}

\keyword{print}
