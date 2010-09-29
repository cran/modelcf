\name{CL_comparts}
\alias{checkParts}
\alias{varInfo}
\alias{countPart}
\title{Comparing partitions (clustering)}

\description{
	\code{checkParts} is an assymetric measure of the matching of P relatively to P_ref.

	The two next indices are symmetric.

	\code{varInfo} computes the variation of information index from Meila article.

	\code{countPart} is a simple counter of matched elements, e.g. the matching level of \code{(1,1,1,2)} and \code{(1,1,2,3)} is 2.
}

\usage{
	checkParts(P, P_ref)

	varInfo(P1, P2)

	countPart(P1, P2)
}

\arguments{
	\item{P,P_ref,P1,P2}{a partition of some data, as outputs by \code{\link{gtclusts}} ; e.g., \code{(1,1,1,1,2,2,2,2,1,1,3,3,3)}}
}

\details{
	All indices are normalized to lie in the range \code{(0,1)}.

	The \code{checkParts} method uses P clusters overlap over P_ref ones to compute an adequation index. It is quite severe, designed for testing of clustering methods.

	The \dQuote{variation of information} index of Meila is a (mathematical) measure between partitions. This is actually a nice property ; see article.
}

\value{
	A real number between 0 and 1, indicating the matching level between the two partitions.
}

\references{
	M. Meila, \bold{Comparing Clusterings}, Statistics Technical Report 418, University of Washington, 2002
}

\examples{
	#comparing the three indices
	P = c(1,1,2,2,2,2,2,2,3,3,3,4,4,4,1,1)
	P_ref = c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
	print(checkParts(P, P_ref))
	print(varInfo(P, P_ref))
	print(countPart(P, P_ref))
}

