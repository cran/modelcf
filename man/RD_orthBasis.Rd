\name{orthBasis}
\alias{linEmb}
\alias{linear_rec}
\title{Around orthonormal bases}

\description{
	\code{linEmb} performs decomposition onto an orthonormal basis amoung functional PCA, wavelets (any filter), Fourier and B-spline basis.

	\code{linear_rec} performs linear reconstruction based on coefficients.
}

\usage{
	linEmb(data, d, linbt="PCA", filt="haar", wvar=TRUE)

	linear_rec(basis, coefs)
}

\arguments{
	\item{data}{matrix of n functions (written as vectors) in rows ; \code{data[i,]} is the i-th D-dimensional function}

	\item{d}{desired number of coefficients ; corresponds to basis resolution, reduced \code{d}-dimensionality}

	\item{linbt}{the type of (linear) orthonormal basis ; \dQuote{PCA} for functional PCA, \dQuote{wav} for wavelets basis, \dQuote{four} for Fourier basis and \dQuote{bsp} for B-spline basis}

	\item{filt}{the desired filter in case of wavelets basis ; choice between EXTREMAL PHASE (daublet): \dQuote{haar}, \dQuote{dX} where X belongs to (4, 6, 8, 10, 12, 14, 16, 18, 20); LEAST ASYMMETRIC (symmlet): \dQuote{sX} where X belongs to (4, 6, 8, 10, 12, 14, 16, 18, 20); BEST LOCALIZED: \dQuote{lX} where X belongs to (2, 4, 6, 14, 18, 20); COIFLET: \dQuote{cX} where X belongs to (6, 12, 18, 24, 30)}

	\item{wvar}{boolean telling if we should select the sub-basis with most variable coefficients}
	
	\item{basis}{orthonormal functions (written as vectors) in rows}
	
	\item{coefs}{matrix of projected coefficients in rows}
}

\value{
	\code{linEmb} returns a list L, with L$embed = matrix of d-dimensional embeddings in rows, L$basis = matrix of orthonormal functions (in rows).

	\code{linear_rec} performs linear reconstruction based on coefficients.
}

\references{
	Functional PCA: J. Ramsay and B. W. Silverman, \bold{Functional Data Analysis}, Springer 2005
	
	Wavelets basis R package used is \emph{wmtsa} available here \url{http://cran.r-project.org/web/packages/wmtsa/index.html}
}

\examples{
	#generate a \dQuote{triginometric} functional dataset
	cosFunc = cos( seq( from=0,to=2*pi,by=2*pi/200 ) )
	sinFunc = sin( seq( from=0,to=2*pi,by=2*pi/200 ) )
	coefs = matrix( runif(200),ncol=2 )
	fdata = coefs \%*\% rbind(cosFunc, sinFunc)
	#plot the two first Fourier functions
	four = linEmb(fdata, 2, "four")
	plotC(four$basis)
	#output the three first PCA functions
	fpca = linEmb(fdata, 3, "PCA")
	plotC(fpca$basis)
}

