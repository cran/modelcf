\name{RD_redDim}
\alias{nlin_redDim}
\alias{nlin_adaptRec}
\title{Dimensionality reduction and associate reconstruction}

\description{
	\code{nlin_redDim} nlin_redDim is a generic method for dimensionality reduction.

	\code{nlin_adaptRec} is a generic method for reconstruction.
	
	For internal use only ; use specific methods directly if you need.
}

\usage{
	nlin_redDim(method, data, d, adn, knn, alpha, knnmin, knnmax,
		tsoft, thlvl, hdth)

	nlin_adaptRec(method, embobj, newEmb)
}

\arguments{
	\item{method}{the dimensionality reduction method (to be) used}

	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{d}{estimated data dimension. It can be estimated using functions \code{\link{dimest1}} or \code{\link{dimest2}}}

	\item{adn}{boolean for adapted point-varying neighborhoods, from Wang et al. article ; in short, the more linear data is around \code{x}, the more \code{x} has neighbors}

	\item{knn}{fixed number of neighbors at each point (used only if \code{adn==FALSE}). If zero, a simple heuristic will determine it around \code{sqrt(nrow(data))}}

	\item{alpha}{fraction of overlapping elements when building the traversal sequence of neighborhoods}

	\item{knnmin}{minimum number of neighbors at each point}

	\item{knnmax}{maximum number of neighbors at each point}
	
	\item{tsoft}{tolerance factor for the visibility graph computation (between 0 and 1 ; should be close to 0)}

	\item{thlvl}{fraction of total elements of data to be embedded using the initial local basis}

	\item{hdth}{\dQuote{hard} threshold, same as above parameter but integer. It defines the maximum level of elements in the Dijkstra graph which will be embedded using the initial local basis. If zero, only \code{thlvl} is considered}

	\item{embobj}{an object as outputs by \code{\link{RML}}, \code{\link{LPcaML}} or \code{\link{GCEM}} functions}

	\item{newEmb}{a new embedding from which the high dimensional object has to be estimated}
}

\value{
	A list with the embedding in \code{$embed}, and some technical parameters for reconstruction.
}

\references{
	T. Lin, H. Zha and S. U. Lee, \bold{Riemannian Manifold Learning for Nonlinear Dimensionality Reduction}, at European Conference on Computer Vision, Graz, Austria 9: 44-55, 2006

	J. Wang, Z. Zhang and H. Zha, \bold{Adaptive Manifold Learning}, in Advances in Neural Information Processing Systems 17: 1473-1480, 2005

	Y. Zhan, J. Yin, G. Zhang and En Zhu, \bold{Incremental Manifold Learning Algorithm Using PCA on Overlapping Local Neighborhoods for Dimensionality Reduction}, at  3rd International Symposium on Advances in Computation and Intelligence 5370: 406-415, 2008
}

