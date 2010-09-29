\name{CL_clustering}
\alias{phclust}
\alias{gtclusts}
\title{Main clustering function}

\description{
	\code{phclust} performs R hierarchical cluster (using hclust()) with Ward linkage, and call cutree() after.

	This function should not be called directly. Use the following one instead.

	\code{gtclusts} main function to cluster data according to any method.
}

\usage{
	phclust(dissims, K)

	gtclusts(method, data, K, d=min(10, ncol(data)), adn=FALSE, knn=0,
		symm=TRUE, weight=FALSE, sigmo=FALSE, alpha=1.0)
}

\arguments{
	\item{method}{the clustering method, to be chosen between \dQuote{HDC} (k-means based on Hitting Times), \dQuote{CTH} (Commute-Time Hierarchic), \dQuote{CTHC} (Commute-Time CHAMELEON), \dQuote{CTKM} (Commute-Time k-means), \dQuote{specH} (\dQuote{spectral-hierarchical} clustering), \dQuote{specKM} (spectral clustering with k-means), \dQuote{CH} (hierarchical clustering), \dQuote{CHC} (CHAMELEON clustering), \dQuote{PCA} (PCA-k-means from Chiou and Li ; see references), \dQuote{KM} (basic k-means)}

	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{dissims}{matrix of dissimilarities (can be simple L2 distances, or more complicated like commute-times)}

	\item{K}{expected number of clusters}

	\item{d}{estimated data dimension (needed only for \dQuote{ACP} method and when \code{adn} is \code{TRUE}). It can be estimated using functions \code{\link{dimest1}} or \code{\link{dimest2}}}

	\item{adn}{boolean for adapted point-varying neighborhoods, from Wang et al. article ; in short, the more linear data is around \code{x}, the more \code{x} has neighbors}

	\item{knn}{fixed number of neighbors at each point ; used only if \code{adn == FALSE}. If zero, a simple heuristic will determine it around \code{sqrt(nrow(data))}}

	\item{symm}{boolean at \code{TRUE} for symmetric similarity matrix (see code. It does not impact much the result}

	\item{weight}{boolean at \code{TRUE} for weighted hitting/commute times, like in the article of Liben-Nowell and Kleinberg}

	\item{sigmo}{boolean at \code{TRUE} for sigmoid commute-time kernel, like in the article of Yen et al.}

	\item{alpha}{parameter controlling the relative importance of clusters' connectivity in CHAMELEON clustering ; usual values range from 0.5 to 2}
}

\details{
	\dQuote{Safe} methods are HDC, CTH, CTKM, CH, PCA and KM. Others could output weird results.

	The spectral clustering is taken from the article of Ng et al., and adapted to work on a possibly disconnected graph

	\code{adn} should not be set when working with small datasets and/or in low dimension (<= 3)

	When \code{sigmo} is set, the sigmoid commute-time kernel (Yen et al.) is computed with a=1. In the paper authors say it need manual tuning.
}

\value{
	An integer vector describing classes (same as \code{kmeans()$cluster} field).
}

\references{
	J-M. Chiou and P-L. Li, \bold{Functional clustering and identifying substructures of longitudinal data}, in Journal of the Royal Statistical Society 69(4): 679-699, 2007

	G. Karypis, E.-H. Han and V. Kumar, \bold{CHAMELEON: A Hierarchical Clustering Algorithm Using Dynamic Modeling}, in IEEE Computer 32(8): 68-75, 1999

	A. Y. Ng, M. Jordan and Y. Weiss, \bold{On Spectral Clustering: Analysis and an algorithm}, at Advances in Neural Information Processing Systems, Vancouver, BC, Canada 14: 849-856, 2002

	D. Liben-Nowell and J. Kleinberg ; \bold{The link-prediction problem for social networks}, in Journal of the American Society for Information Science and Technology 58(7): 1019-1031, 2007

	J. Wang, Z. Zhang and H. Zha, \bold{Adaptive Manifold Learning}, in Advances in Neural Information Processing Systems 17: 1473-1480, 2005

	L. Yen, D. Vanvyve, F. Wouters, F. Fouss, M. Verleysen and M. Saerens, \bold{Clustering using a random-walk based distance measure}, at Symposium on Artificial Neural Networks 13: 317-324, Bruges, Belgium, 2005

	L. Yen, F. Fouss, C. Decaestecker, P. Francq and M. Saerens, \bold{Graph nodes clustering with the sigmoid commute-time kernel: A comparative study}, in Data \& Knowledge Engineering 68(3): 338-361, 2009
}

\examples{
	#generate a mixture of three gaussian data sets
	data = rbind( matrix(rnorm(200,mean=2,sd=0.5),ncol=2),
	    matrix(rnorm(200,mean=4,sd=0.5),ncol=2),
	    matrix(rnorm(200,mean=6,sd=0.5),ncol=2) )
	#cluster it using k-means
	km = gtclusts("KM", data, 3)
	#and using Commute-Time Hierarchic clustering
	ct = gtclusts("CTH", data, 3, knn=20, symm=FALSE)
	#plot results
	plotPts(data, cl=km)
	plotPts(data, cl=ct)
}

