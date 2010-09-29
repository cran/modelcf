\name{Z_predict}
\alias{predictClassif}
\alias{predictRegress}
\alias{predict.modelcf}
\title{Predictions for some models}

\description{
	\code{predictClassif} estimates the label of an object x.

	\code{predictRegress} estimates the output y for an input x.

	These two last functions should not be used directly. Prefer calling specific methods from some R package.

	\code{predict.modelcf} estimates the output curve y for an input vector x, using a model built by the \code{\link{fmetam}} function.
}

\usage{
	predictRegress(model,newIn_s)

	predictClassif(model,newIns)

	\method{predict}{modelcf}(object, x, verb = FALSE, \dots)
}

\arguments{
	\item{model}{a classification or regression model, as output by \code{\link{learnClassif}} or \code{\link{learnRegress}}}

	\item{newIn_s, newIns}{a matrix of (testing) input vectors in rows}

	\item{object}{a \pkg{modelcf} model, output of \code{\link{fmetam}}}

	\item{x}{a matrix of n input vectors in rows, which can be given as a R matrix or a text file. x[i,] is the i-th p-dimensional input.}

	\item{verb}{TRUE for printing what is going on. A further release will allow to choose levels of verbosity}

	\item{\ldots}{unused (for compatibility with generic method \code{predict})}
}

\value{
	\code{predictClassif} (resp. \code{predictRegress}) returns a vector of integer (resp. real) values.

	\code{predict.modelcf} return a matrix of curves in rows, one for each testing example.
}

\examples{
	#get the first artificial dataset and build a standard model of it
	#using 250 training samples
	data(datacf)
	trainInds = sample(1:300, 250)
	m = fmetam(dataIn[trainInds,],dataOut1[trainInds,],d=3,wcl=FALSE,mdim="linear")
	#get the predicted curves
	pred = predict.modelcf(m, dataIn[-trainInds,])
	#plot the (L1) error between real and predicted curves
	plotC(dataOut1[-trainInds,] - pred)
}

