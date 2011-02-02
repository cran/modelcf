\name{redDim}
\alias{nlin_redDim}
\alias{nlin_adaptRec}
\title{Dimensionality reduction and associate reconstruction}

\description{
	\code{nlin_redDim} nlin_redDim is a generic method for dimensionality reduction.

	\code{nlin_adaptRec} is a generic method for reconstruction.
	
	For internal use only ; use specific methods directly if you need.
}

\usage{
	nlin_redDim(method, data, d, adn, k, alpha, trcv,
		kmin, kmax, tsoft, thlvl, hdth, advhr)

	nlin_adaptRec(method, embobj, newEmb)
}

\arguments{
	\item{method}{the dimensionality reduction method (to be) used, to be chosen between \code{\link{RML}} and \code{\link{LPcaML}}}

	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{d}{estimated data dimension (e.g. through functions from \code{\link{dimension}} file)}

	\item{adn}{string for adapted point-varying neighborhoods. "none" for no adaptivity, "adbas" for simple local PCA based neighborhoods (see code), "ad1" for the Zhan et al. method, and "ad2" for Wang et al. method. In short, the more linear data is around \code{x}, the more \code{x} has neighbors}

	\item{k}{fixed number of neighbors at each point (used only if \code{adn=="none"}). If zero, a simple heuristic will determine it around \code{sqrt(nrow(data))}}

	\item{alpha}{fraction of overlapping elements when building the traversal sequence of neighborhoods}

	\item{trcv}{fraction of total examples on which a model is trained during cross-validation procedures}

	\item{kmin}{minimum number of neighbors at each point}

	\item{kmax}{maximum number of neighbors at each point}
	
	\item{tsoft}{tolerance factor for the visibility graph computation (between 0 and 1 ; should be close to 0)}

	\item{thlvl}{fraction of total elements of data to be embedded using the initial local basis}

	\item{hdth}{\dQuote{hard} threshold, same as above parameter but integer. It defines the maximum level of elements in the Dijkstra graph which will be embedded using the initial local basis. If zero, only \code{thlvl} is considered}

	\item{advhr}{if TRUE, the heuristic for RML's last step is based on estimated graph distances from minimal connectivity graph ; if FALSE, the heuristic use euclidian distances}

	\item{embobj}{an object as output by \code{\link{RML}} or \code{\link{LPcaML}} functions}

	\item{newEmb}{a new embedding from which the high dimensional object has to be estimated}
}

\value{
	A list with the embedding in \code{$embed}, and some technical parameters for reconstruction.
}

