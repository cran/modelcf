\name{RML}
\alias{RML}
\alias{RML_rec}
\title{Riemannian Manifold Learning}

\description{
	\code{RML} embeds data in the \code{d}-dimensional space using the Riemannian Manifold Learning method from the Lin et al. article.

	\code{RML_rec} inverses the above procedure, reconstructing a curve (or any high dimensional vector) from its low-dimensional representation.
}

\usage{
	RML(data, d, kmin=0, kmax=0, tsoft=0.1,
		thlvl=0.3, hdth=2, advhr=FALSE)

	RML_rec(RLout, newEmb)
}

\arguments{
	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{d}{estimated data dimension (e.g. through functions from \code{\link{dimension}} file)}

	\item{kmin}{minimum number of neighbors at each point}

	\item{kmax}{maximum number of neighbors at each point}

	\item{tsoft}{tolerance factor for the visibility graph computation (between 0 and 1 ; should be close to 0)}

	\item{thlvl}{fraction of total elements of data to be embedded using the initial local basis}

	\item{hdth}{\dQuote{hard} threshold, same as above parameter but integer. It defines the maximum level of elements in the Dijkstra graph which will be embedded using the initial local basis. If zero, only \code{thlvl} is considered}

	\item{advhr}{if TRUE, the heuristic for RML's last step is based on estimated graph distances from minimal connectivity graph ; if FALSE, the heuristic use euclidian distances}

	\item{RLout}{an object as output by \code{RML} function}

	\item{newEmb}{a new embedding from which the high dimensional object has to be estimated}
}

\details{
	The algorithm works in two main steps :

	\enumerate{
	\item{An origin vector \code{y0} is determined, and its neighbors are embedded by projection onto a local tangent basis.}
	\item{For further away elements \code{y}, we first find the predecessor \code{yp} of \code{y} on a shortest path from \code{y0}, and the \code{yp} neighbors written \code{yi1,...,yik}. The core idea then is to preserve (as much as possible) angles \code{y-yp-yij} to get the embedding \code{z}.}
	}

	The reconstruction function \code{RML_rec} does exactly the same thing but from low-dimensional space to high-dimensional one. For better explanations, see the article.
}

\value{
	A list with the embedding in \code{$embed}, and some technical parameters for reconstruction.
}

\references{
	T. Lin, H. Zha and S. U. Lee, \bold{Riemannian Manifold Learning for Nonlinear Dimensionality Reduction}, at European Conference on Computer Vision, Graz, Austria 9: 44-55, 2006

	T. Lin and H. Zha, \bold{Riemannian Manifold Learning}, in IEEE Transactions on Pattern Analysis and Machine Intelligence 30 (5): 796-809, 2008
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
	emb = RML(sw, 2, kmin=15, kmax=30)$embed
	plotPts(emb, cl=cols)
}

