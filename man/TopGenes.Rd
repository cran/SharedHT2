\name{TopGenes}
\alias{TopGenes}
\title{Gene lists linked to Gene Cards}
\description{
  \code{TopGenes} creates a genelist sorted on values of a chosen statistic.
}
\usage{
  TopGenes(obj, by = "EB", FDR = 0.05, allsig = FALSE, n.g = 20, 
           browse = FALSE, search.url = genecards, subset.contrasts = NULL,
           subset.cols.sig = NULL, on.alpha = NULL, off.alpha = NULL,
           fold.change = NULL, out.values = c("mean", "tstat", "both"),
           path = "", file = "") 
}
\arguments{
  \item{obj}{An object of class \code{fit.n.data} returned by \code{EB.Anova}}
  \item{by}{Specify "EB" or "naive".  Do you want to use the empirical bayes 
    or the naive variant of the statistic. Defaults to "EB".}
  \item{FDR}{The false discovery fate (FDR) to use in the Benjamini-Hochberg
    (BH) stepdown procedure.}
  \item{allsig}{Set to \code{TRUE} if you only want to see the genes 
    meeting the BH criterion. Defaults to FALSE}
  \item{n.g}{If \code{allsig} is \code{FALSE}, you must specify the number
    of top genes that you want to see.}
  \item{browse}{set to \code{TRUE} if you want the results displayed in the HTML
    browser with gene identifiers linked to the GeneCards database at the 
    Weizmann institute. Defaults to FALSE. See examples below}


  \item{search.url}{should contain a url href to search an online database
    when a gene identifier is appended onto the end. Defaults to 'genecards',
    included in this package, for searching the GeneCards database at the
    Weizmann institute.}
  \item{subset.cols.sig, subset.contrasts, on.alpha, off.alpha, fold.change}{
    If desired, the list obtained from the BH FDR stepdown procedure can be
    filtered further based upon nested T tests and/or fold change applied to
    a supplied list of filtering contrasts.  The arguments \code{subset.contrasts}
    and \code{subset.cols.sig} are used to control this post-BHFDR filtering.
    The argument \code{subset.contrasts} is a matrix with column dimension equal
    to the number of group labels and arbitrary row dimension. The user may test
    contrasts of the original experimental groups by placing one contrast in each
    row. The argument \code{subset.cols.sig} is a vector of length equal to the row
    dimension of the contrasts matrix. The components of this vector are given the
    values 1 or -1 depending upon whether filtering should be based upon the contrast
    being 'on' or 'off'.  Contrasts are combined via 'OR' for 'on' and 'AND' for 'off'
    so that the filtered list consists of all genes for which all 'off' contrasts are
    'off' and one or more of the 'on' contrasts is 'on'.  The conditions 'on'/'off' are
    determined via nested T-statistics and or fold change.  To use nested T-filtering,
    specify a type I error in the argument \code{on.alpha} and if you are using any
    'off' contrasts, a type II error in the argument \code{off.alpha}. To use
    fold change filtering specify a fold change threshold in the argument
    \code{fold.change}}
  \item{out.values}{Use this to specify what information you want displayed in the
    output table, by making a selection from  c("mean", "tstat", "both").  The default is
    "mean"}
  \item{file}{If \code{browse} is set to \code{TRUE} and if you want to 
    save the html file then supply a filename, \code{file=foobgwb.html}}
  \item{path}{If specifying a \code{file} name above, optionally you
    may specify a \code{path} name. In the unix implementation, specifying
    a file name without a path writes to the current working directory.
    This feature is not yet supported in Windows. Use \code{path} argument
    to explicitly specify a directory}
}

\value{
  A list containing two components
  \item{table}{sorted genelist having the following columns:
    \code{RowNum}, the row number from the original unsorted data frame. Useful when you
    happen to know that the first 100 genes are true positives and the rest are not
    (as is the case with the supplied dataset, 'SimAffyDat'),  \code{GeneId},
    which is taken from the row component of dimnames in the original dataframe.  As such,
    it is very usefull to name your rows using the affy gene identifiers as these are
    searchable in the GeneCards database, \code{'NAME1'} ...\code{'NAMEd'}, the group means
    corresponding to the group of each name,  \code{'TYPE'.stat}, the statistic used in
    the sort. This is determined by user choices at two levels. First, when computations
    are performed inside the call to 'EB.Anova', the specification of 'Var.Struct', which
    defaults to "general" with alternate value "simple" determins whether the
    multivariate test, i.e. the Hotelling T-squared, (HT2) or the univariate F test (UT2)
    are computed. In both cases, both the empirical Bayes and standard variants are
    computed. Next after the computation is completed, when the user requests the sorted
    genelist, the specification of the arguement 'by', either "EB" (default) or "naive"
    determines which of the two computed statistics is used to perform the sort.  Thus
    'TYPE' assumes one of four values 'ShHT2', 'HT2', 'ShUT2', or 'UT2', where HT2 versus
    UT2 depend upon the specification of 'Var.Struct' used in 'EB.Anova' (general versus
    naive), while the prefix 'Sh' refers to shrinkgage variance, which depends upon how
    sorting is to be done, specified in the argument \code{by} ('EB' for shrinkage and
    'naive' otherwise). The next column, \code{'TYPE'.p-val}, is the corresponding
    p-value under the corresponding model.  The final column, \code{FDR.stepdown='FDR'},
    contains the BH criterion values which are computed as 'FDR * rank/Ngenes'}
  \item{call}{The function call producing the result.}
}



\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{EB.Anova}}, \code{\link{EBfit}}, \code{\link{SimAffyDat}},
         \code{\link{TopGenes}}, \code{\link{SimNorm.IG}},
         \code{\link{SimMVN.IW}}, \code{\link{SimMVN.mxIW}}, 
         \code{\link{SimOneNorm.IG}}, \code{\link{SimOneMVN.IW}},
         \code{\link{SimOneMVN.mxIW}}
}

\note{The print method for class \code{fit.n.data} contains a call to 
  \code{TopGenes} so that calls to \code{TopGenes} and to \code{print}
  work almost the same, except the latter produces a model fit summary
  as well.}

\examples{
#  NOTE:  Full browser functionality includes (i) links to the NCI 'GENECARDS'
#  database and (ii) gene name floats above the affy identifier when mouse
#  pointer is over the link.  Right now this is specific to Affymetrix arrays
#  more functionality may be added later.  The hot-links feature, (i), works
#  provided the dataset is given affy identifiers for rownames.  These are
#  the default if you used the Bioconductor package 'affy' suite of preprocessing
#  routines.  The floating gene names on mouse-over function, (ii), works provided
#  you attach the annotation name as an attribute to your dataset.
#  For exammple, suppose you used the 'affy' suite to do pre-processing
#  and have obtained the exprSet, 'MyEset', and will use  MYDAT <- exprs(MyEset)
#  in calls to EB.Anova. Then function (i) is automatic, and to get function (ii) to
#  work requires the following simple command:
#
#  attr(MYDAT, "annotation") <- annotation(MyEset)
#
#
#
# The included example dataset is a simulated hgu95av2 Affymetrix oligonucleotide
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
