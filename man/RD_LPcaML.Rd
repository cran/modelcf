\name{LPcaML}
\alias{LPcaML}
\alias{LPcaML_rec}
\title{Local PCA Manifold Learning}

\description{
	\code{LPcaML} embeds data in the \code{d}-dimensional space using the Local PCA Manifold Learning method from the Zhan et al. article.

	\code{LPcaML_rec} inverses the above procedure, reconstructing a curve (or any high dimensional vector) from its low-dimensional representation.
}

\usage{
	LPcaML(data, d, adn="none", k=0, alpha=0.5, trcv=0.7)

	LPcaML_rec(LPout, newEmb)
}

\arguments{
	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{d}{estimated data dimension (e.g. through functions from \code{\link{dimension}} file)}

	\item{adn}{string for adapted point-varying neighborhoods. "none" for no adaptivity, "adbas" for simple local PCA based neighborhoods (see code), "ad1" for the Zhan et al. method, and "ad2" for Wang et al. method. In short, the more linear data is around \code{x}, the more \code{x} has neighbors}

	\item{k}{fixed number of neighbors at each point (used only if \code{adn==FALSE}). If zero, a simple heuristic will determine it around \code{sqrt(nrow(data))}}

	\item{alpha}{fraction of overlapping elements when building the traversal sequence of neighborhoods}

	\item{trcv}{fraction of total examples on which a model is trained during cross-validation procedures}
	
	\item{LPout}{an object as output by \code{LPcaML} function}

	\item{newEmb}{a new embedding from which the high dimensional object has to be estimated}
}

\details{
	The algorithm works in two main steps :

	\enumerate{
	\item{A traversal sequence of (overlapping) local neighborhoods is constructed, and a simple PCA is computed in each neighborhood.}
	\item{The reduced coordinates are then computed step by step, by optimizing an affine transformation matrix on the overlap between two neighborhoods.}
	}

	The reconstruction function \code{LPcaML_rec} first find the right neighborhood, then apply inverse affine transformation. For better explanations, see the article.
}

\value{
	A list with the embedding in \code{$embed}, and some technical parameters for reconstruction.
}

\references{
	Y. Zhan, J. Yin, G. Zhang and E. Zhu, \bold{Incremental Manifold Learning Algorithm Using PCA on Overlapping Local Neighborhoods for Dimensionality Reduction}, at  3rd International Symposium on Advances in Computation and Intelligence 5370: 406-415, 2008
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
	emb = LPcaML(sw, 2, alpha=0.7)$embed
	plotPts(emb, cl=cols)
}

