\name{paste}
\alias{paste}
\alias{\%,\%}
\title{The paste operator}
\description{A binary operator shortcut for paste(x,y)}
\usage{
x \%,\% y
}
\arguments{
  \item{x}{a character string}
  \item{y}{a character string}
}
 
\value{
  paste(x, y, sep="")
}
\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\examples{
  
  "log2.grp" \%,\% (1:2)
  c(outer("log2.grp" \%,\% (1:2), ".n" \%,\% (1:3), FUN="\%,\%"))

}

\keyword{utilities}
