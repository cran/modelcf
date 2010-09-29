\name{CR_planExp}
\alias{xtr_plan1}
\alias{xtr_plan2}
\title{Extract training set}

\description{
	\code{xtr_plan1} and \code{xtr_plan2} extract a training set from n samples.
}

\usage{
	xtr_plan1(data, knn, trcv)

	xtr_plan2(data, knn, trcv)
}

\arguments{
	\item{data}{matrix of \code{n} vectors in rows. \code{data[i,]} is the i-th m-dimensional vector}

	\item{knn}{fixed number of neighbors at each point to build the training set in cross-validation procedure}

	\item{trcv}{fraction of total examples on which a model is trained during cross-validation procedure.}
}

\details{
	\code{xtr_plan1} searches for the maximum local variance points.
	\code{xtr_plan2} searches for the maximum local density points.
}

\value{
	An integer vector giving the indices of the selected plan.
}

\examples{
	#generate a mixture of three gaussian data set
	data = rbind( matrix(rnorm(200,mean=2,sd=0.5),ncol=2),
	    matrix(rnorm(200,mean=4,sd=0.5),ncol=2),
	    matrix(rnorm(200,mean=6,sd=0.5),ncol=2) )
	#extract a plan of type 1
	x1 = xtr_plan1(data, 20, 0.7)
	cols = rep(1,300) ; cols[x1] = 2
	plotPts(data,cl=cols)
	#extract a plan of type 2
	x2 = xtr_plan2(data, 20, 0.7)
	cols = rep(1,300) ; cols[x2] = 2
	plotPts(data,cl=cols)
}

