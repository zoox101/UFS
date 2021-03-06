% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UFS.R
\name{variance_explained}
\alias{variance_explained}
\title{Variance Explained}
\usage{
variance_explained(df, idxs, show = TRUE, ...)
}
\arguments{
\item{df}{The dataframe to plot.}

\item{idxs}{Specific indices to search from the dataframe.}

\item{show}{Shows the plot. Defaults to TRUE.}
}
\description{
Calculates the variance explained by a set of vectors.
}
\examples{
variance_explained(senate, c(1,2,3), main="Variance Explained By [1,2,3]")
variance_explained(t(senate), c(1,2,3), main="Variance Explained By [1,2,3]")
}
\keyword{explained}
\keyword{variance,}
