\name{regression}
\alias{learnRegress}
\alias{optimParams_regress}
\title{Statistical learning (regression)}

\description{
	\code{learnRegress} builds a regression object (see code for details).

	\code{optimParams_regress} optimize parameters for the chosen method.

	These two methods should not be called directly. Using the specific technique inside its own package is a better idea.
}

\usage{
	learnRegress(x, y, method, params, stred)

	optimParams_regress(x, y, method, k, d, trcv, verb)
}

\arguments{
	\item{x}{matrix of \code{n} input vectors in rows. \code{x[i,]} is the i-th p-dimensional input}

	\item{y}{matrix of \code{n} outputs in rows. \code{y[i,]} is the i-th m-dimensional output}

	\item{method}{regression method, to be chosen between \dQuote{PPR} (Projection Pursuit Regression), \dQuote{rforest} (random forests), \dQuote{kNN, fkNN} (Nadaraya-Watson, after dimensionality reduction or not), \dQuote{lPCA} (local PCA regression, without dimensionality reduction), \dQuote{GP} (gaussian processes), \dQuote{SVR} (Support Vector Regression)}

	\item{params}{vector of parameters for the chosen method}

	\item{stred}{boolean at \code{TRUE} for standardize outputs \code{y}}

	\item{k}{fixed number of neighbors at each point to build the training set in cross-validation procedure}

	\item{d}{estimated outputs dimensionality ; relevant only if data has not been reduced}

	\item{trcv}{fraction of total examples on which a model is trained during cross-validation procedure.}

	\item{verb}{\code{TRUE} for printing what is going on.}
}

\value{
	\code{learnRegress} returns a regression object (internal specifications).

	\code{optimParams_regress} returns a vector of optimized parameters for the chosen method.
}

