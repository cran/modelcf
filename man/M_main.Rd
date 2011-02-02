\name{modeling}
\alias{fmetam_1cl}
\alias{fmetam}
\alias{fm_resids}
\alias{nfoldcv}
\title{Build and validate the functional outputs model}

\description{
	\code{fmetam_1cl} is a subroutine to do the dimensionality reduction step. Internal use only.

	\code{fmetam} is the main method to build a model, using clustering and dimensionality reduction.

	\code{fm_resids} will call \code{fmetam} twice, first on basic data and then on residuals for a better fit.

	\code{nfoldcv} builds and tests several models with fixed parameters and ramdomly generated training sets (cross-validation).
}

\usage{
	fmetam_1cl(x, y, d, mdim, adnRD, knnRD, linbt, filt, wvar, alpha,
	kmin, kmax, tsoft, thlvl, hdth, advhr, mreg, ppts, stred, trcv, verb)

	fmetam(x, y, d=0, mclust="CTH", mclass="kNN", redy=FALSE, adnCC="none",
	knnCC=0, wcl=TRUE, iclusts=rep(0,nrow(y)), symm=TRUE, weight=FALSE, 
	sigmo=FALSE, minszcl=30, maxcl=Inf, taus=0.8, Ns=10, tauc=0.8, Nc=10,
	mdim="linear", adnRD="none", knnRD=0, linbt="PCA", filt="haar",
	wvar=TRUE, alpha=0.5, kmin=0, kmax=0, tsoft=0.1, thlvl=0.3, hdth=2, 
	advhr=FALSE, mreg="PPR", ppts=FALSE, stred=TRUE, trcv = 0.7, verb = TRUE)

	fm_resids(x, y, d=0, mclust="CTH", mclass="kNN", redy=FALSE, adnCC="none",
	knnCC=0, wcl=TRUE, iclusts=rep(0,nrow(y)), symm=TRUE, weight=FALSE, 
	sigmo=FALSE, minszcl=30, maxcl=Inf, taus=0.8, Ns=10, tauc=0.8, Nc=10,
	mdim1="linear", mdim2="RML", adnRD="none", knnRD=0, linbt="PCA", 
	filt="haar", wvar=TRUE, alpha=0.5, kmin=0, kmax=0, tsoft=0.1, thlvl=0.3, 
	hdth=2, advhr=FALSE, mreg1="PPR", mreg2="PPR", ppts=FALSE, stred=TRUE, 
	trcv = 0.7, verb = TRUE)

	nfoldcv(x, y, d=0, single=TRUE, mclust="CTH", mclass="kNN", redy=FALSE,
	adnCC="none", knnCC=0, wcl=TRUE, symm=TRUE, weight=FALSE, sigmo=FALSE,
	minszcl=30, maxcl=Inf, taus=0.8, Ns=10, tauc=0.8, Nc=10,
	mdim1="linear", mdim2="RML", adnRD="none", knnRD=0, linbt="PCA", filt="haar",
	wvar=TRUE, alpha=0.5, kmin=0, kmax=0, tsoft=0.1, thlvl=0.3, hdth=2,
	advhr=FALSE, mreg1="PPR", mreg2="PPR", ppts=FALSE, stred=TRUE, trcv = 0.7,
	loo = FALSE, nfold=100, nhold=10, verb = TRUE, plotc=TRUE)
}

\arguments{
	\item{x}{matrix of \code{n} input vectors in rows, given as a R matrix or filename. \code{x[i,]} is the i-th p-dimensional input}

	\item{y}{matrix of \code{n} discretized outputs in rows, given as a R matrix or filename. \code{y[i,]} is the i-th D-dimensional output}

	\item{d}{estimated (real) outputs dimensionality (should be far less than \code{D})}

	\item{single}{boolean telling if the model will be composite (base + residuals)}

	\item{mdim, mdim1, mdim2}{the dimensionality reduction method (to be) used for base model (1) or residuals (2) : choice between \dQuote{linear} for orthonormal basis, \dQuote{RML} for Riemannian Manifold Learning, \dQuote{LPcaML} for Local PCA Manifold Learning and \dQuote{GCEM} for Global Coordination of Exponential Maps}

	\item{adnRD}{string for adapted point-varying neighborhoods in dimensionality reduction. "none" for no adaptivity, "adbas" for simple local PCA based neighborhoods (see code), "ad1" for the Zhan et al. method, and "ad2" for Wang et al. method. In short, the more linear data is around \code{x}, the more \code{x} has neighbors}

	\item{knnRD}{fixed number of neighbors at each point for dimensionality reduction (used only if \code{adnRD=="none"}). If zero, a simple heuristic will determine it around \code{sqrt(nrow(data))}. Irrelevant if \code{mdim=="RML"}}

	\item{linbt}{the type of (linear) orthonormal basis ; \dQuote{PCA} for functional PCA, \dQuote{wav} for wavelets basis, \dQuote{four} for Fourier basis and \dQuote{bsp} for B-spline basis}

	\item{filt}{the desired filter in case of wavelets basis ; choice between EXTREMAL PHASE (daublet): \dQuote{haar}, \dQuote{dX} where X belongs to (4, 6, 8, 10, 12, 14, 16, 18, 20); LEAST ASYMMETRIC (symmlet): \dQuote{sX} where X belongs to (4, 6, 8, 10, 12, 14, 16, 18, 20); BEST LOCALIZED: \dQuote{lX} where X belongs to (2, 4, 6, 14, 18, 20); COIFLET: \dQuote{cX} where X belongs to (6, 12, 18, 24, 30)}

	\item{wvar}{boolean telling if we should select the sub-basis with most variable coefficients}

	\item{alpha}{fraction of overlapping elements when building the traversal sequence of neighborhoods}

	\item{kmin}{minimum number of neighbors at each point}

	\item{kmax}{maximum number of neighbors at each point}

	\item{tsoft}{tolerance factor for the visibility graph computation (between 0 and 1 ; should be close to 0)}

	\item{thlvl}{fraction of total elements of data to be embedded using the initial local basis}

	\item{hdth}{\dQuote{hard} threshold, same as above parameter but integer. It defines the maximum level of elements in the Dijkstra graph which will be embedded using the initial local basis. If zero, only \code{thlvl} is considered}

	\item{advhr}{if TRUE, the heuristic for RML's last step is based on estimated graph distances from minimal connectivity graph ; if FALSE, the heuristic use euclidian distances}

	\item{mreg, mreg1, mreg2}{regression method to use (1 for base model, 2 for residuals) ; choice between between \dQuote{PPR} (Projection Pursuit Regression), \dQuote{rforest} (random forests), \dQuote{kNN, fkNN} (Nadaraya-Watson, after dimensionality reduction or not), \dQuote{lPCA} (local PCA regression, without dimensionality reduction), \dQuote{GP} (gaussian processes), \dQuote{SVR} (Support Vector Regression)}

	\item{ppts}{TRUE for pointwise regression}
	
	\item{stred}{TRUE for standardized outputs}
	
	\item{trcv}{fraction of total examples on which a model is trained during cross-validation procedures}
	
	\item{mclust}{the clustering method, to be chosen between \dQuote{HDC} (k-means based on Hitting Times), \dQuote{CTH} (Commute-Time Hierarchic), \dQuote{CTKM} (Commute-Time k-means), \dQuote{spec} (spectral clustering), \dQuote{CH} (hierarchical clustering), \dQuote{PCA} (PCA-k-means from Chiou and Li ; see references), \dQuote{KM} (basic k-means)}

	\item{mclass}{type of classifier to use in the prediction accuracy step ; choice between \dQuote{kNN} (k-nearest-neighbors), \dQuote{ctree} (classification tree), \dQuote{RDA} (Regularized Discriminant Analysis), \dQuote{rforest} (random forests), \dQuote{SVM} (Support Vector Machines). Only the first two were intensively tested}

	\item{redy}{boolean telling if the outputs should be reduced (with PCA) as a preprocessing step}

	\item{adnCC}{same as \code{adnRD} above, for clustering}

	\item{knnCC}{fixed number of neighbors at each point in clustering ; used only if \code{adnCL == FALSE}}

	\item{wcl}{FALSE for disable clustering step ; can be useful for comparison purposes}
	
	\item{iclusts}{imposed clustering, like (1,1,2,2,2,1) (if known by user ; used in the \code{fm_resids} method)}
	
	\item{symm}{boolean at \code{TRUE} for symmetric similarity matrix (see code. It does not impact much the result}

	\item{weight}{boolean at \code{TRUE} for weighted hitting/commute times, like in the article of Liben-Nowell and Kleinberg}

	\item{sigmo}{boolean at \code{TRUE} for sigmoid commute-time kernel, like in the article of Yen et al.}

	\item{minszcl}{minimum size for a cluster. This is interesting to not allow too small clusters for the regression stage ; recommanded values are above 30-50}

	\item{maxcl}{maximum number of clusters ; \code{Inf} stands for \dQuote{no limit}, i.e. determined by stability-prediction loops only}

	\item{taus}{threshold for stability check ; value between 0 (every method accepted) and 1 (only ultra-stable method accepted). Recommanded between 0.6 and 0.9}

	\item{Ns}{number of stability runs before averaging results (the higher the better, although slower..)}

	\item{tauc}{threshold for prediction accuracy check (after subsampling) ; value between 0 (every clustering accepted) and 1 (only \dQuote{well separated} clusters accepted). Recommanded between 0.6 and 0.9}

	\item{Nc}{number of partitions predictions runs before averaging results (same remark as for \code{Ns} above)}

	\item{loo}{TRUE for leave-one-out cross-validation} 

	\item{nfold}{number of cross-validation loops to run}
	
	\item{nhold}{number of curves to hold in the training step for cross-validation}

	\item{verb}{\code{TRUE} for printing what is going on. A further release will allow to choose levels of verbosity}
	
	\item{plotc}{TRUE for plotting current Q2 curves at each step}
}

\details{
	If code{d} argument is left unspecified (0), it will be estimated using Farahmand et al. algorithm.

	The algorithm in \code{fmetam} works in three main steps :

	\enumerate{
	\item{Optional clustering of intputs-outputs.}
	\item{Dimensionality reduction in each outputs cluster.}
	\item{Statistical learning "inputs --> reduced coordinates".}
	}

	The \code{\link{predict.modelcf}} function then computes the associated reconstruction "recuced coordinates --> curves".
}

\value{
	\code{fmetam_1cl} and \code{fmetam} return a list of relevant parameters for internal use.
	
	\code{nfoldcv} returns a list with the following attributes:
	\itemize{
	\item{curves = predicted curves (only in leave-one-out mode);}
	\item{MSE = (average) mean squares error curve for the model chosen;}
	\item{stMSE = corresponding standard deviation;}
	\item{pvar = (average) mean squares error curve for the training mean model;}
	\item{stvar = corresponding standard deviation;}
	\item{Q2 = Q2 error curve (should be above 0 and close to 1);}
	\item{stQ2 = corresponding standard deviation;}
	\item{ssclust = measure of clusters' sizes homogeneity (>=0, should be as small as possible);}
	\item{snclust = histogram vector of the number of clusters found over the runs; e.g., \code{(0,0,32,78,0,...,0)} means 78 runs with 4 clusters and 32 runs with 3 clusters.}
	}
	NOTE: standard deviations cannot be accurate if \code{nfold} parameter is too small. Value around 100 or above is recommended.
}

\examples{
	data(datacf)
	#plot curves of the dataset
	plotC(dataOut1)
	plotC(dataOut2)

	#build a standard model of the first dataset using 250 training samples
	trainInds = sample(1:300, 250)
	m = fmetam(dataIn[trainInds,],dataOut1[trainInds,],d=3,wcl=FALSE,mdim="linear")
	# print the model
	print(m)
	#get the predicted curves
	pred = predict.modelcf(m, dataIn[-trainInds,])
	#get and plot error estimators
	errs = fperrors(pred,dataOut1[-trainInds,],colMeans(dataOut1[trainInds,]))
	plot(errs$MSE, type="l", ylab="MSE")
	plot(1-errs$MSE/errs$pvar, type="l", ylim=c(0,1), ylab="Q2")

	# run cross validation for the second dataset
	\dontrun{nf = nfoldcv(dataIn,dataOut2,d=3,wcl=FALSE,mdim="linear",plotc=FALSE)}
	nf = nfoldcv(dataIn[1:200,],dataOut2[1:200,],d=3,wcl=FALSE,mdim1="linear",
    nfold=10,plotc=FALSE) #for speed
	# plot MSE +/- standard deviation
	rg = range(nf$MSE-nf$stMSE,nf$MSE+nf$stMSE)
	plot(nf$MSE-nf$stMSE,type="l",lwd=3,col=4,ylim=rg); par(new=TRUE)
	plot(nf$MSE+nf$stMSE,type="l",lwd=3,col=4,ylim=rg); par(new=TRUE)
	plot(nf$MSE,type="l",lwd=3,col=1,ylim=rg)
	# plot Q2 +/- standard deviation
	rg = c(-0.5, 1.5)
	plot(nf$Q2-nf$stQ2,type="l",lwd=3,col=4,ylim=rg); par(new=TRUE)
	plot(nf$Q2+nf$stQ2,type="l",lwd=3,col=4,ylim=rg); par(new=TRUE)
	plot(nf$Q2,type="l",lwd=3,ylim=rg)
}

