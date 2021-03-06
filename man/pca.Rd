% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UFS.R
\name{pca}
\alias{pca}
\title{Principle Component Analysis}
\usage{
pca(df, n, show = TRUE, ...)
}
\arguments{
\item{df}{The dataframe to plot.}

\item{n}{The number of dimensions to keep.}

\item{show}{Displays the plot. Defaults to TRUE.}
}
\value{
\item{transformation_matrix}{The transformation matrix.}

\item{var_explained}{The explained variance by each component.}
}
\description{
Runs PCA and plots the variance explained by a set of column vectors.
}
\examples{
pca(senate, 2, main="Bill PCA Variance Explained")
pca(t(senate), 2, main="Senator PCA Variance Explained")
}
\keyword{pca}
