\name{CL_kmeans}
\alias{km_PCA}
\alias{km_dists}
\title{k-means like functions}

\description{
	\code{km_dists} = k-means based on a distance matrix.

	\code{km_PCA} = generalization of classical k-means for functional case, by Chiou and Li.
}

\usage{
	km_dists(distm, K, nstart=10, maxiter=100)

	km_PCA(data, K, d, simplif=TRUE, maxiter=50)
}

\arguments{
	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{distm}{matrix of distances (can be simple L2 distances, or more complicated like commute-times)}

	\item{K}{expected number of clusters}

	\item{d}{estimated data dimension. It can be estimated using functions \code{\link{dimest1}} or \code{\link{dimest2}}}

	\item{simplif}{boolean at \code{TRUE} for simplified algorithm, without leave-one-out SVD's (actually very costly)}

	\item{nstart}{number of algorithm runs with random initialization}

	\item{maxiter}{maximum number of iterations within one algorithm run}
}

\details{
	The k-means using a distances matrix is exactly the same algorithm as classical k-means, except for the choice of centroids, which must belong to the dataset.

	The PCA-k-means algorithm replaces the centroids by centroids \bold{plus} local basis functions obtained by (functional) PCA. The closeness to a cluster is computed relatively to this full system, instead of a centroid only. Apart from this point, the algorithm is similar to k-means ; but more general. The \code{simplif} argument allows or not a simplification avoiding very costly leave-one-out procedure, (re)computing local basis after slight data change. It can be switched off without fears for big enough clusters (say, more then a few dozens).
}

\value{
	An integer vector describing classes (same as \code{kmeans()$cluster} field).
}

\references{
	J-M. Chiou and P-L. Li, \bold{Functional clustering and identifying substructures of longitudinal data}, in Journal of the Royal Statistical Society 69(4): 679-699, 2007
}

\examples{
	#generate a mixture of three gaussian data set, and compute distances
	data = rbind( matrix(rnorm(200,mean=2,sd=0.5),ncol=2),
	    matrix(rnorm(200,mean=4,sd=0.5),ncol=2),
	    matrix(rnorm(200,mean=6,sd=0.5),ncol=2) )
	dists = as.matrix(dist(data))
	#cluster using k-means
	km = km_dists(dists, 3)
	#plot result
	plotPts(data, cl=km)
	#and using km_PCA clustering after artificial functional transformation
	sPoints = seq(from=0,to=2*pi,by=2*pi/200)
	cosFunc = cos(sPoints)
	sinFunc = sin(sPoints)
	fdata = as.matrix(data[,1]) \%*\% cosFunc + as.matrix(data[,2]^2) \%*\% sinFunc
	kp = km_PCA(fdata, 3, 2)
	#plot result
	plotC(fdata, cl=kp)
}

