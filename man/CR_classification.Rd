\name{CR_classification}
\alias{learnClassif}
\alias{optimParams_classif}
\title{Building a classifier}

\description{
	\code{learnClassif} builds a classifier object (see code for details).

	\code{optimParams_classif} optimize parameters for the chosen method.

	These two methods should not be called directly. Using the specific technique inside its own package is a better idea.
}

\usage{
	learnClassif(x, y, method, params)

	optimParams_classif(x, y, method, knn, trcv)
}

\arguments{
	\item{x}{matrix of \code{n} input vectors in rows. \code{x[i,]} is the i-th p-dimensional input}

	\item{y}{matrix of \code{n} outputs in rows. \code{y[i,]} is the i-th m-dimensional output}

	\item{method}{classification method, to be chosen between \dQuote{kNN} (k-nearest-neighbors), \dQuote{ctree} (classification trees), \dQuote{RDA} (Regularized Discriminant Analysis), \dQuote{rforest} (random forests), \dQuote{SVM} (Support Vector Machines)}

	\item{params}{vector of parameters for the chosen method}

	\item{knn}{fixed number of neighbors at each point to build the training set in cross-validation procedure}

	\item{trcv}{fraction of total examples on which a model is trained during cross-validation procedure.}
}

\value{
	\code{learnClassif} returns a classifier object (internal specifications).

	\code{optimParams_classif} returns a vector of optimized parameters for the chosen method.
}

