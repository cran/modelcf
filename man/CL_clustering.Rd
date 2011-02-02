\name{clustering}
\alias{clustering}
\alias{phclust}
\alias{gtclusts}
\alias{findK_gtclusts}
\alias{gtclusts_inout}
\title{Main clustering functions}

\description{
	\code{phclust} performs R hierarchical cluster (using hclust()) with Ward linkage, and call cutree() after.

	This function should not be called directly. Use the following one instead.

	\code{gtclusts} main function to cluster data according to any method.

	\code{findK_gtclusts} is a procedure to determine the number of classes (and associate partitioning).

	\code{gtclusts_inout} calls the previous method one on outputs, and then on each inputs cluster (main procedure).
}

\usage{
	phclust(dissims, K)

	gtclusts(method, data, K, d=min(5,ncol(data)), adn="none", knn=0,
		symm=TRUE, weight=FALSE, sigmo=FALSE)

	findK_gtclusts(x, y, method, d=min(10, ncol(x)), adn="none", knn=0,
		symm=TRUE, weight=FALSE, sigmo=FALSE, minszcl=30,
		maxcl=Inf, mclass="kNN", taus=0.8, Ns=10, tauc=0.8, Nc=10,
		trcv=0.7, nstagn=10)

	gtclusts_inout(x, y, method, d=min(10, ncol(x)), redy=FALSE, adn="none",
		knn=0, symm=TRUE, weight=FALSE, sigmo=FALSE, minszcl=30,
		maxcl=Inf, mclass="kNN", taus=0.8, Ns=10, tauc=0.8, Nc=10,
		trcv=0.7, verb=TRUE, nstagn=10)
}

\arguments{
	\item{method}{the clustering method, to be chosen between \dQuote{HDC} (k-means based on Hitting Times), \dQuote{CTH} (Commute-Time Hierarchic), \dQuote{CTKM} (Commute-Time k-means), \dQuote{spec} (spectral clustering), \dQuote{CH} (hierarchical clustering), \dQuote{PCA} (PCA-k-means from Chiou and Li ; see references), \dQuote{KM} (basic k-means)}

	\item{data}{matrix of \code{n} vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}
	
	\item{x}{matrix of \code{n} input vectors in rows. \code{x[i,]} is the i-th p-dimensional input}

	\item{y}{matrix of \code{n} discretized outputs in rows. \code{y[i,]} is the i-th D-dimensional output}

	\item{dissims}{matrix of dissimilarities (can be simple L2 distances, or more complicated like commute-times)}

	\item{K}{expected number of clusters}

	\item{d}{estimated (real) outputs dimensionality (should be far less than \code{D}) ; useful only if one of the following parameters is set: \code{redy},\code{adn},\code{method=="ACP"}. It can be estimated using functions from \code{\link{dimension}} file}

	\item{adn}{string for adapted point-varying neighborhoods. "none" for no adaptivity, "adbas" for simple local PCA based neighborhoods (see code), "ad1" for the Zhan et al. method, and "ad2" for Wang et al. method. In short, the more linear data is around \code{x}, the more \code{x} has neighbors}

	\item{knn}{fixed number of neighbors at each point ; used only if \code{adn=="none"}. If zero, a simple heuristic will determine it around \code{sqrt(nrow(data))}}

	\item{symm}{boolean at \code{TRUE} for symmetric similarity matrix (see code. It does not impact much the result}

	\item{weight}{boolean at \code{TRUE} for weighted hitting/commute times, like in the article of Liben-Nowell and Kleinberg}

	\item{sigmo}{boolean at \code{TRUE} for sigmoid commute-time kernel, like in the article of Yen et al.}

	\item{redy}{boolean telling if the outputs should be reduced (with PCA) as a preprocessing step}

	\item{minszcl}{minimum size for a cluster. This is interesting to not allow too small clusters for the regression stage ; recommanded values are above 30-50}

	\item{maxcl}{maximum number of clusters ; \code{Inf} stands for \dQuote{no limit}, i.e. determined by stability-prediction loops only}

	\item{mclass}{type of classifier to use in the prediction accuracy step ; choice between \dQuote{kNN} (k-nearest-neighbors), \dQuote{ctree} (classification tree), \dQuote{RDA} (Regularized Discriminant Analysis), \dQuote{rforest} (random forests), \dQuote{SVM} (Support Vector Machines). Only the first two were intensively tested}

	\item{taus}{threshold for stability check ; value between 0 (every method accepted) and 1 (only ultra-stable method accepted). Recommanded between 0.6 and 0.9}

	\item{Ns}{number of stability runs before averaging results (the higher the better, although slower..)}

	\item{tauc}{threshold for prediction accuracy check (after subsampling) ; value between 0 (every clustering accepted) and 1 (only \dQuote{well separated} clusters accepted). Recommanded between 0.6 and 0.9}

	\item{Nc}{number of partitions predictions runs before averaging results (same remark as for \code{Ns} above)}

	\item{trcv}{fraction of total examples on which a model is trained during cross-validation procedures.}

	\item{verb}{\code{TRUE} for printing what is going on. A further release will allow to choose levels of verbosity.}

	\item{nstagn}{number of allowed stages (increasing the number of clusters K) without added clusters (if \code{minszcl} is large enough small clusters may end being merged).}
}

\details{
	\code{adn} should not be set when working with small datasets and/or in low dimension (<= 3).

	When \code{sigmo} is set, the sigmoid commute-time kernel (Yen et al.) is computed with a=1. In the paper authors say it need manual tuning.

	The algorithm for simultaneous estimate of \code{K} and clustering works in two main steps :

	\enumerate{
	\item{subsample original data in data1 and data2, then cluster both, and measure similarity between partitions at the intersection using the variation of information index of Meila article.}
	\item{subsample a training set Tr in \code{[1,n]} where n is the number of data rows, then subsample a set S which must contain \code{[1,n] \ Tr}. Cluster both sets, and use Tr to predict labels of the testing set. Finally compare the partitions using simple \dQuote{matching counter} after renumbering (with the hungarian algorithm).}
	}

	Both are repeated \code{Ns, Nc} times to get accurate estimators.
	We stop when these estimators fall below the thresholds \code{taus, tauc}, and return corresponding partition.
}

\value{
	An integer vector describing classes (same as \code{kmeans()$cluster} field).
}

\references{
	J-M. Chiou and P-L. Li, \bold{Functional clustering and identifying substructures of longitudinal data}, in Journal of the Royal Statistical Society 69(4): 679-699, 2007

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
	ct = gtclusts("CTH", data, 3, k=20, symm=FALSE)
	#plot results
	plotPts(data, cl=km)
	plotPts(data, cl=ct)
	
	#generate a (smaller) mixture of three gaussian data sets
	inData = rbind( matrix(rnorm(60,mean=2,sd=0.5),ncol=2),
	    matrix(rnorm(60,mean=4,sd=0.5),ncol=2),
	    matrix(rnorm(60,mean=6,sd=0.5),ncol=2) )
	#build artificial corresponding outputs
	sPoints = seq(from=0,to=2*pi,by=2*pi/200)
	cosFunc = cos(sPoints)
	sinFunc = sin(sPoints)
	outData = as.matrix(inData[,1]) \%*\% cosFunc + as.matrix(inData[,2]^2) \%*\% sinFunc
	#partition inputs-outputs using Commute-Time Hierarchic clustering
	ct = gtclusts_inout(inData, outData, "CTH", k=20, minszcl=20, mclass="kNN",
	    taus=0.7, Ns=10, tauc=0.7, Nc=10)
	#plot results, inputs then outputs
	plotPts(inData, cl=ct)
	plotC(outData, cl=ct)
}

