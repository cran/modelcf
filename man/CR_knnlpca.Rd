\name{CR_knnlpca}
\alias{knnPredict}
\alias{lpcaPredict}
\title{Prediction for no-dimensionality-reduction-methods}

\description{
	\code{knnPredict and lpcaPredict} predict the response(s) for the given input(s), using \code{fkNN} and \code{lPCA} models.

	These functions should not be called directly. Use \code{\link{predictRegress}} instead.
}

\usage{
	knnPredict(x_train, y_train, x, knn)

	lpcaPredict(x_train, y_train, x, knn, stzouts=TRUE)
}

\arguments{
	\item{x_train}{matrix of \code{n} training input vectors in rows. \code{x_train[i,]} is the i-th p-dimensional training input}

	\item{y_train}{matrix of \code{n} training outputs in rows. \code{y_train[i,]} is the i-th m-dimensional output}

	\item{x}{matrix (or vector) of \code{q} testing input vectors in rows. \code{x[i,]} is the i-th p-dimensional testing input}

	\item{knn}{fixed number of neighbors at each point}

	\item{stzouts}{boolean at \code{TRUE} for standardize outputs \code{y_train} (after internal dimensionality reduction)}
}

\value{
	A matrix of predictions, in same format as \code{y_train}.
}

