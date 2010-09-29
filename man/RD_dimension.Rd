\name{RD_dimension}
\alias{locdim1}
\alias{dimest1}
\alias{dimest2}
\title{Dimension estimation}

\description{
	\code{locdim1} estimates the local dimension at some point (a row of the distance matrix).

	\code{dimest1} estimates the intrinsic dimension of \code{data}, following the algorithm of Farahmand et al.

	\code{dimest2} uses the \code{\link{RML}} graph to estimate dimension through its simplices as indicated by Lin et al. (this is really slow).
}

\usage{
	locdim1(x, knn)

	dimest1(data, knn)

	dimest2(data, knnmin, knnmax, tsoft=0.1)
}

\arguments{
	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{x}{a row of the distance matrix computed from \code{data}}

	\item{knn}{fixed number of neighbors at each point}

	\item{knnmin}{minimum number of neighbors at each point}

	\item{knnmax}{maximum number of neighbors at each point}

	\item{tsoft}{tolerance factor for the visibility graph computation (between 0 and 1 ; should be close to 0)}
}

\value{
	An integer equals to the estimated dimension.
}

\references{
	A. M. Farahmand, C. Szepesvari and J-Y. Audibert, \bold{Manifold-adaptive dimension estimation}, at 24th International Conference on Machine Learning 227: 265-272, 2007

	T. Lin, H. Zha and S. U. Lee, \bold{Riemannian Manifold Learning for Nonlinear Dimensionality Reduction}, at European Conference on Computer Vision, Graz, Austria 9: 44-55, 2006
}

\examples{
	#generate a swissroll dataset
	n = 300 ; h = 3
	phi = runif(n, min=0, max=2*pi)
	z = runif(n, min=0, max=h)
	sw = cbind( phi*cos(phi), phi*sin(phi), z )

	#estimate dimension
	print(dimest1(sw, 20))
	\dontrun{print(dimest2(sw, 15, 25)) #WARNING: very very slow!}
}

