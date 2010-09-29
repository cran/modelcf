\name{CL_findK_Clust}
\alias{findK_gtclusts}
\alias{gtclusts_inout}
\title{Clustering input-output data}

\description{
	\code{findK_gtclusts} is a procedure to determine the number of classes (and associate partitioning).

	\code{gtclusts_inout} calls the previous method one on outputs, and then on each inputs cluster (main procedure).
}

\usage{
	findK_gtclusts(x, y, method, d=min(10, ncol(data)), adn=FALSE, knn=0,
		symm=TRUE, weight=FALSE, sigmo=FALSE, alpha=1.0, minszcl=30,
		maxcl=Inf, mclass="kNN", taus=0.95, Ns=10, tauc=0.95, Nc=10,
		trcv=0.7, nstagn=10)

	gtclusts_inout(x, y, method, d=min(10, ncol(data)), redy=TRUE, adn=FALSE,
		knn=0, symm=TRUE, weight=FALSE, sigmo=FALSE, alpha=1.0, minszcl=30,
		maxcl=Inf, mclass="kNN", taus=0.95, Ns=10, tauc=0.95, Nc=10,
		trcv=0.7, verb=TRUE, nstagn=10)
}

\arguments{
	\item{x}{matrix of \code{n} input vectors in rows. \code{x[i,]} is the i-th p-dimensional input}

	\item{y}{matrix of \code{n} discretized outputs in rows. \code{y[i,]} is the i-th D-dimensional output}

	\item{method}{clustering method, to be chosen between \dQuote{HDC} (k-means based on Hitting Times), \dQuote{CTH} (Commute-Time Hierarchic), \dQuote{CTHC} (Commute-Time CHAMELEON), \dQuote{CTKM} (Commute-Time k-means), \dQuote{specH} (\dQuote{spectral-hierarchical} clustering), \dQuote{specKM} (spectral clustering with k-means), \dQuote{CH} (hierarchical clustering), \dQuote{CHC} (CHAMELEON clustering), \dQuote{PCA} (ACP-k-means from Chiou and Li ; see references), \dQuote{KM} (basic k-means)}

	\item{d}{estimated (real) outputs dimensionality (should be far less than \code{D}) ; useful only if one of the following parameters is set: \code{redy},\code{adn},\code{method=="ACP"}}

	\item{redy}{boolean telling if the outputs should be reduced (with PCA) as a preprocessing step}

	\item{adn}{boolean for adapted point-varying neighborhoods, from Wang et al. article ; in short, the more linear data is around \code{x}, the more \code{x} has neighbors}

	\item{knn}{fixed number of neighbors at each point ; used only if \code{adn == FALSE}}

	\item{symm}{boolean at \code{TRUE} for symmetric similarity matrix (see code. It does not impact much the result}

	\item{weight}{boolean at \code{TRUE} for weighted hitting/commute times, like in the article of Liben-Nowell and Kleinberg}

	\item{sigmo}{boolean at \code{TRUE} for sigmoid commute-time kernel, like in the article of Yen et al.}

	\item{alpha}{parameter controlling the relative importance of clusters' connectivity in CHAMELEON clustering ; usual values range from 0.5 to 2}

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
	The algorithm works in two main steps :

	\enumerate{
	\item{subsample original data in data1 and data2, then cluster both, and measure similarity between partitions at the intersection using the variation of information index of Meila article.}
	\item{subsample a training set Tr in \code{[1,n]} where n is the number of data rows, then subsample a set S which must contain \code{[1,n] \ Tr}. Cluster both sets, and use Tr to predict labels of the testing set. Finally compare the partitions using simple \dQuote{matching counter} after renumbering (with the hungarian algorithm).}
	}

	Both are repeated \code{Ns, Nc} times to get accurate estimators.
	We stop when these estimators fall below the thresholds \code{taus, tauc}, and return corresponding partition.
}

\value{
	An integer vector describing classes (same as \code{kmeans()$cluster} field). The number of clusters is equal to its maximum value.
}

\references{
	J-M. Chiou and P-L. Li, \bold{Functional clustering and identifying substructures of longitudinal data}, in Journal of the Royal Statistical Society 69(4): 679-699, 2007

	G. Karypis, E.-H. Han and V. Kumar, \bold{CHAMELEON: A Hierarchical Clustering Algorithm Using Dynamic Modeling}, in IEEE Computer 32(8): 68-75, 1999

	A. Y. Ng, M. Jordan and Y. Weiss, \bold{On Spectral Clustering: Analysis and an algorithm}, at Advances in Neural Information Processing Systems, Vancouver, BC, Canada 14: 849-856, 2002

	D. Liben-Nowell and J. Kleinberg ; \bold{The link-prediction problem for social networks}, in Journal of the American Society for Information Science and Technology 58(7): 1019-1031, 2007

	M. Meila, \bold{Comparing Clusterings}, Statistics Technical Report 418, University of Washington, 2002

	J. Wang, Z. Zhang and H. Zha, \bold{Adaptive Manifold Learning}, in Advances in Neural Information Processing Systems 17: 1473-1480, 2005

	L. Yen, D. Vanvyve, F. Wouters, F. Fouss, M. Verleysen and M. Saerens, \bold{Clustering using a random-walk based distance measure}, at Symposium on Artificial Neural Networks 13: 317-324, Bruges, Belgium, 2005

	L. Yen, F. Fouss, C. Decaestecker, P. Francq and M. Saerens, \bold{Graph nodes clustering with the sigmoid commute-time kernel: A comparative study}, in Data \& Knowledge Engineering 68(3): 338-361, 2009
}

\examples{
	#generate a mixture of three gaussian data sets
	\dontrun{inData = rbind( matrix(rnorm(200,mean=2,sd=0.5),ncol=2),
	    matrix(rnorm(200,mean=4,sd=0.5),ncol=2),
	    matrix(rnorm(200,mean=6,sd=0.5),ncol=2) ) #use this data matrix}
	inData = rbind( matrix(rnorm(60,mean=2,sd=0.5),ncol=2),
	    matrix(rnorm(60,mean=4,sd=0.5),ncol=2),
	    matrix(rnorm(60,mean=6,sd=0.5),ncol=2) ) #for speed
	#build artificial corresponding outputs
	sPoints = seq(from=0,to=2*pi,by=2*pi/200)
	cosFunc = cos(sPoints)
	sinFunc = sin(sPoints)
	outData = as.matrix(inData[,1]) \%*\% cosFunc + as.matrix(inData[,2]^2) \%*\% sinFunc
	#partition only outputs using hierarchical clustering
	ch = findK_gtclusts(inData, outData, "CH", knn=20, minszcl=20, mclass="kNN",
	    taus=0.7, Ns=10, tauc=0.7, Nc=10)
	#plot result
	plotC(outData, cl=ch)
	#partition inputs-outputs using Commute-Time Hierarchic clustering
	ct = gtclusts_inout(inData, outData, "CTH", knn=20, minszcl=20, mclass="kNN",
	    taus=0.7, Ns=10, tauc=0.7, Nc=10)
	#plot results, inputs then outputs
	plotPts(inData, cl=ct)
	plotC(outData, cl=ct)
}

