% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UFS.R
\name{scree_plot}
\alias{scree_plot}
\title{Scree Plot}
\usage{
scree_plot(xs, ys = NULL, cumulate = FALSE)
}
\arguments{
\item{xs}{The first series to plot.}

\item{ys}{An additional second series to plot.}

\item{cumulate}{Plots a cumulative scree plot. Defaults to false.}
}
\description{
Plots a scree plot of one or two vectors
}
\examples{
scree_plot(pca(senate, 10, FALSE), gvm(senate, 10, FALSE), TRUE)
scree_plot(pca(t(senate), 10, FALSE), gvm(t(senate), 10, FALSE), TRUE)
}
\keyword{plot}
\keyword{scree,}
