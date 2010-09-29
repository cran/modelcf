\name{CL_chameleon}
\alias{updateDissims}
\alias{chameleon}
\title{CHAMELEON clustering}

\description{
	\code{updateDissims} auxiliar function to update dissimilarities between one cluster and every other.

	\code{chameleon} main function to cluster data according to CHAMELEON method, described in the article of Karypis et al. given in reference.

	These two functions should not be called directly. Use \code{\link{gtclusts}} instead.
}

\usage{
	updateDissims(data, dissims, clusts, flags, index, alpha)

	chameleon(data, dissims, K, alpha)
}

\arguments{
	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{dissims}{matrix of dissimilarities (can be simple L2 distances, or more complicated like commute-times)}

	\item{K}{expected number of clusters}

	\item{alpha}{parameter controlling the relative importance of clusters' connectivity ; usual values range from 0.5 to 2}

	\item{clusts}{the vector of current clusters (something like \code{(1 1 1 2 2 3 1 3 3 4 4)} ... describing classes)}

	\item{flags}{boolean vector ; \code{flags[i] == TRUE} if \code{i} is the representative element for its cluster}

	\item{index}{index from which dissimilarities must be computed}
}

\value{
	An integer vector describing classes (same as \code{kmeans()$cluster} field).
}

\references{
	George Karypis, Eui-Hong Han and Vipin Kumar, CHAMELEON: A Hierarchical Clustering Algorithm Using Dynamic Modeling, in IEEE Computer 32(8): 68-75, 1999
}

