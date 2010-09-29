\name{Z_mixpred}
\alias{getcoefc}
\alias{mixpredf}
\title{Mixing functional models}

\description{
	Functions to define a mixture of already created models.

	\code{getcoefc} returns a curve matching the maximums given by user (to facilitate models mixing).

	\code{mixpredf} takes several models as arguments, and mix them after calling \code{\link{predict.modelcf}}. This allows to benefit from different kinds of models.
}

\usage{
	getcoefc(D, inds, maxs=c(), rgs=c())

	mixpredf(mods, coefs, x, verb = FALSE)
}

\arguments{
	\item{D}{outputs dimensionality (usually a few hundreds)}

	\item{inds}{(strictly) positive integer vector of desired local maximums locations}

	\item{maxs}{positive real vector of desired local maximums amplitudes}

	\item{rgs}{minimum number of neighbors at each point}

	\item{mods}{a list of \pkg{modelcf} models, outputs of \code{\link{fmetam}}}

	\item{coefs}{a list of curves (same length as training outputs), which are taken as mixture coefficients (see details below)}

	\item{x}{matrix of n testing input vectors in rows ; \code{x[i,]} is the i-th m-dimensional testing input vector}

	\item{verb}{\code{TRUE} for printing what is going on. A further release will allow to choose levels ov verbosity.}
}

\details{
	\code{getcoefc} outputs a piecewise constant function, which locally constant parts are centered around indices given in \code{inds}. The (integer) width of each locally constant part is given by the \code{rgs} vector argument ; if not provided, the width is taken constant, equals to the maximum value which avoid overlapping. \code{maxs} indicates the amplitude of each local maximum (piecewise constant), and will equals \code{(1,1,1,1,\ldots)} if not provided.
}

\value{
	\code{getcoefc} returns a sampled curve (with \code{D} values).

	\code{mixpredf} returns a model prediction (matrix with curves in rows) ; same output format as \code{\link{predict.modelcf}}.
}

\examples{
	#get the first artificial dataset and build three different models of it
	#using 250 training samples
	data(datacf)
	trainInds = sample(1:300, 250)
	m1 = fmetam(dataIn[trainInds,], dataOut1[trainInds,], d=3, wcl=FALSE,
	    mdim="linear")
	m2 = fmetam(dataIn[trainInds,], dataOut1[trainInds,], d=3, wcl=FALSE,
	    mdim="RML", knnmin=15, knnmax=25)
	m3 = fmetam(dataIn[trainInds,], dataOut1[trainInds,], d=3, wcl=FALSE,
	    mreg="fkNN")
	#mix the three, giving \dQuote{first third} weight to the first,
	#\dQuote{second third} weight to the second
	#and \dQuote{third third} weight to the third one
	mix = mixpredf(list(m1,m2,m3), list(c(rep(1,66),rep(0,134)),
	    c(rep(0,66),rep(1,67),rep(0,67)),c(rep(0,133),rep(1,67))),
	    dataIn[-trainInds,], verb=TRUE)
	#plot the (L1) error between real and predicted curves
	plotC(dataOut1[-trainInds,] - mix)
}

