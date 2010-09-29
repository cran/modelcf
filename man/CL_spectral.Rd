\name{CL_spectral}
\alias{spec_emb}
\alias{spec_clust}
\title{Spectral clustering}

\description{
	\code{spec_emb} embeds data into \code{R^d}, using symmetric laplacian (see Ng et al. article).

	\code{spec_clust} clusters data into \code{K} classes using spectral embedding.

	These two functions should not be called directly. Use \code{\link{gtclusts}} instead.
}

\usage{
	spec_emb(data, K, d, adn, knn)

	spec_clust(method, data, K, d, adn, knn)
}

\arguments{
	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{method}{the clustering (sub-)method, to be chosen between \dQuote{specH} (\dQuote{spectral-hierarchical} clustering) and \dQuote{specKM} (spectral clustering with k-means)}

	\item{K}{expected number of clusters}

	\item{d}{estimated data dimension, e.g. as output by functions \code{\link{dimest1}} or \code{\link{dimest2}}}

	\item{adn}{boolean for adapted point-varying neighborhoods, from Wang et al. article ; in short, the more linear data is around \code{x}, the more \code{x} has neighbors}

	\item{knn}{fixed number of neighbors at each point ; used only if \code{adn == FALSE}}
}

\details{
	The \code{spec_clust} procedures works also for disconnected graphs by merging connected components (if more than \code{K}), or subdivising largest clusters when \code{K} is bigger than the number of connected components.
}

\value{
	\code{spec_emb} returns a matrix embedding the data in rows.

	\code{spec_clust} returns an integer vector describing classes (same as \code{kmeans()$cluster} field).
}

\references{
	A. Y. Ng, M. Jordan and Y. Weiss, \bold{On Spectral Clustering: Analysis and an algorithm}, at Advances in Neural Information Processing Systems, Vancouver, BC, Canada 14: 849-856, 2002

	J. Wang, Z. Zhang and H. Zha, \bold{Adaptive Manifold Learning}, in Advances in Neural Information Processing Systems 17: 1473-1480, 2005
}

