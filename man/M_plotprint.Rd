\name{M_printPlot}
\alias{plotC}
\alias{plotPts}
\alias{print.modelcf}
\title{Printing and plotting utilitiy functions}

\description{
	\code{plotC} plots a matrix of curves (in rows).

	\code{plotPts} plots a set of 2D points given by the column numbers in a matrix.

	\code{print.modelcf} prints some relevant parameters of a constructed model (as output by \code{\link{fmetam}}).
}

\usage{
	plotC(data, cl=rep(1,nrow(data)), rg=c(min(data),max(data)), \ldots)

	plotPts(data, cols=c(1,2), cl=rep(1,nrow(data)), \ldots)

	\method{print}{modelcf}(x, \ldots)
}

\arguments{
	\item{data}{matrix of n vectors (\dQuote{or functions}) in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{cl}{an integer vector with R colors to be applied to each row}

	\item{rg}{the range on y axis in case of functions drawing}

	\item{cols}{the two selected columns in case of points plotting}

	\item{x}{a model as output by \code{\link{fmetam}}}

	\item{\ldots}{any other relevant graphical parameter(s)}
}

\examples{
	#plot first artificial dataset
	data(datacf)
	plotC(dataOut1)
	#generate a mixture of three gaussian data set
	data = rbind( matrix(rnorm(200,mean=2,sd=0.5),ncol=2),
	    matrix(rnorm(200,mean=4,sd=0.5),ncol=2),
	    matrix(rnorm(200,mean=6,sd=0.5),ncol=2) )
	#cluster it using k-means
	km = gtclusts("KM", data, 3)
	#plot result
	plotPts(data, cl=km)
}

