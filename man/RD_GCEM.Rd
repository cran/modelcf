\name{RD_GCEM}
\alias{GCEM}
\alias{GCEM_rec}
\title{Global Coordination of Exponential Maps}

\description{
	\code{GCEM} embeds data in the \code{d}-dimensional space using a mixing of the Local PCA Manifold Learning method from the Zhan et al. article, and the Riemannian Manifold Learning method from the Lin et al. article.

	\code{GCEM_rec} inverses the above procedure, reconstructing a curve (or any high dimensional vector) from its low-dimensional representation.
}

\usage{
	GCEM(data, d, adn=FALSE, knn=0, alpha=0.5, tsoft=0.1,
		thlvl=0.3, hdth=0)

	GCEM_rec(GCout, newEmb)
}

\arguments{
	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{d}{estimated data dimension. It can be estimated using functions \code{\link{dimest1}} or \code{\link{dimest2}}}

	\item{adn}{boolean for adapted point-varying neighborhoods, from Wang et al. article ; in short, the more linear data is around \code{x}, the more \code{x} has neighbors}

	\item{knn}{fixed number of neighbors at each point (used only if \code{adn==FALSE}). If zero, a simple heuristic will determine it around \code{sqrt(nrow(data))}}

	\item{alpha}{fraction of overlapping elements when building the traversal sequence of neighborhoods}

	\item{tsoft}{tolerance factor for the visibility graph computation (between 0 and 1 ; should be close to 0)}

	\item{thlvl}{fraction of total elements of data to be embedded using the initial local basis}

	\item{hdth}{\dQuote{hard} threshold, same as above parameter but integer. It defines the maximum level of elements in the Dijkstra graph which will be embedded using the initial local basis. If zero, only \code{thlvl} is considered}

	\item{GCout}{an object as output by \code{GCEM} function}

	\item{newEmb}{a new embedding from which the high dimensional object has to be estimated}
}

\details{
	The algorithm works exactly as \code{LPcaML} algorithm, but the local PCA coordinates are replaced by \dQuote{local} RML coordinates. It is very experimental, and currently does not work as well as expected.
}

\value{
	A list with the embedding in \code{$embed}, and some technical parameters for reconstruction.
}

\references{
	T. Lin, H. Zha and S. U. Lee, \bold{Riemannian Manifold Learning for Nonlinear Dimensionality Reduction}, at European Conference on Computer Vision, Graz, Austria 9: 44-55, 2006

	J. Wang, Z. Zhang and H. Zha, \bold{Adaptive Manifold Learning}, in Advances in Neural Information Processing Systems 17: 1473-1480, 2005

	Y. Zhan, J. Yin, G. Zhang and En Zhu, \bold{Incremental Manifold Learning Algorithm Using PCA on Overlapping Local Neighborhoods for Dimensionality Reduction}, at  3rd International Symposium on Advances in Computation and Intelligence 5370: 406-415, 2008
}

\examples{
	#generate a swissroll dataset
	n = 300 ; h = 3
	phi = runif(n, min=0, max=2*pi)
	z = runif(n, min=0, max=h)
	#::set colors
	rSize = 64
	r = rainbow(rSize)
	cols = r[pmin(floor((rSize/(2.0*pi))*phi)+1,rSize)]
	#end set colors::
	sw = cbind( phi*cos(phi), phi*sin(phi), z )

	#launch algorithm and visualize result
	emb = GCEM(sw, 2, alpha=0.7)$embed
	plotPts(emb, cl=cols)
}

