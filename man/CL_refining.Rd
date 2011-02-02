\name{refining}
\alias{reordering}
\alias{fusion_smcl}
\alias{mergeToK}
\title{Rearrangement of clusters}

\description{
	\code{reordering} changes the clusters numerotation to use all the integers from 1 to K.

	\code{fusion_smcl} merges clusters until no one has size inferior than \code{minszcl} argument.

	\code{mergeToK} merges clusters given through its arguments until there are exactly K classes.
}

\usage{
	reordering(clusts)

	fusion_smcl(data, clusts, minszcl)

	mergeToK(data, clusts, K)
}

\arguments{
	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{clusts}{a partition of the data, as outputs by \code{\link{gtclusts}} ; e.g., \code{(1,1,1,1,2,2,2,2,1,1,3,3,3)}}

	\item{K}{expected number of clusters}

	\item{minszcl}{minimum size of a cluster}
}

\value{
	An integer vector describing classes (same as \code{kmeans()$cluster} field).
}

\examples{
	#on an artificial dataset
	data = matrix(runif(300),ncol=3)
	clusts = gtclusts("KM",data,10)
	print(clusts)
	#fusion clusters of size >=20
	print(reordering(fusion_smcl(data,clusts,20)))
	#merge until 3 clusters
	print(mergeToK(data,clusts,3))
}

