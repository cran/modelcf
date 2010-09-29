\name{M_errors}
\alias{fperrors}
\title{Empirical error estimators}

\description{
	\code{fperrors} estimates the error of a model on a specific testing set. It computes MSE errors indicators, by comparing predictions to real curves.
}

\usage{
	fperrors(ypred, yreal, mntrain)
}

\arguments{
	\item{ypred}{matrix of the predicted functions in rows (D sample points / columns)}

	\item{yreal}{matrix of the expected functions (same format as ypred above)}

	\item{mntrain}{mean curve of training outputs}
}

\value{
	A list with MSE values for the model, and the constant estimator (equals to the training mean). The corresponding attributes are named respectively \dQuote{MSE} and \dQuote{pvar}.
}

\examples{
	#get the first artificial dataset and build a standard model of it
	#using 250 training samples
	data(datacf)
	trainInds = sample(1:300, 250)
	m = fmetam(dataIn[trainInds,],dataOut1[trainInds,],d=3,wcl=FALSE,mdim="linear")
	#get the predicted curves and errors
	pred = predict.modelcf(m, dataIn[-trainInds,])
	errs = fperrors(pred,dataOut1[-trainInds,],colMeans(dataOut1[trainInds,]))
	#plot the MSE and Q2 error curves
	plot(errs$MSE, type="l", ylab="MSE")
	plot(1-errs$MSE/errs$pvar, type="l", ylim=c(0,1), ylab="Q2")
}

