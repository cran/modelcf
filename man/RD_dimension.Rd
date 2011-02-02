\name{dimension}
\alias{dimension}
\alias{dest_PCA}
\alias{dest_clust}
\alias{dest_fara}
\alias{dest_pett}
\alias{dest_levi}
\alias{dest_unbl}
\alias{dest_RML}
\alias{dest_rgrl}
\alias{dest_sliv}
\title{Dimension estimation}

\description{
    \code{dest_PCA} and \code{dest_clust} use the PCA locally to estimate dimension, and then average the results. The second one determine local regions using k-means clustering (Bruske and Sommer, 1998), while the first can be considered as a fast suboptimal version of the algorithm by Fan et al. (2010).

	\code{dest_pett}, \code{dest_fara}, \code{dest_levi} and \code{dest_unbl} estimate the intrinsic dimension of \code{data}, following respectively the algorithm of Pettis et al. (1979), Farahmand et al. (2007), Levina and Bickel (2005) and a debiased version of this last one, by MacKay and Ghahramani. All these methods are based on a relation between the dimension and density of data.

    \code{dest_RML} and \code{dest_rgrl} implement an idea from the papers of Lin et al. (2006), using \dQuote{non flat} simplices to evaluate dimension. The second one is a regularization attempt, which has not proven effective yet. \code{dest_sliv} implements a variation on an idea by Cheng and Chiu (2009), simplified although heavier to run.
}

\usage{
    dest_PCA(data, k, thvar=0.01)

    dest_clust(data, nclusts, thvar=0.01)

    dest_pett(data, kmax)

    dest_levi(data, k)

    dest_unbl(data, k)

    dest_fara(data, k)

    dest_RML(data, kmax, N=10, tsoft=0.0)

    dest_rgrl(data, kmax, N=10, alpha=3)
	
	dest_sliv(data, k, N=10000, thtest=0.05)
}

\arguments{
	\item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

	\item{k}{fixed number of neighbors at each point}

	\item{thvar}{expected threshold on explained variance (between 0 and 1 ; should be close to 0)}

    \item{nclusts}{number of cells to be obtained by the k-means algorithm}

    \item{kmax}{maximum number of neighbors at each point}

    \item{N}{number of Monte-Carlo loops}

    \item{tsoft}{tolerance factor for the visibility graph computation (between 0 and 1 ; should be close to 0)}

    \item{alpha}{regularization parameter ; a positive value which give (decreasing) weighted emphasize to increasing values of \code{tsoft}. When \code{alpha} becomes large, the weights associated with \dQuote{high} values of \code{tsoft} fade rapidly}
	
	\item{thtest}{threshold on the p-value of the statistical test for densities adequation}
}

\value{
	An integer equals to the estimated dimension.
}

\references{
    J. Bruske and G. Sommer, \bold{Intrinsic dimension estimation with optimally topology preserving maps}, in IEEE Transactions on Pattern Analysis and Machine Intelligence 20: 572-575, 1998

	M. Fan, N. Gu, H. Qiao and B. Zhang, \bold{Intrinsic dimension estimation of data by principal component analysis}, submitted for publication, 2010.

    K. W. Pettis, T. A. Bailey, A. K. Jain and R. C. Dubes, \bold{An Intrinsic Dimensionality Estimator from Near-Neighbor Information}, in IEEE Transactions on Pattern Analysis and Machine Intelligence 1 (1): 25-37, 1979

	A. M. Farahmand, C. Szepesvari and J-Y. Audibert, \bold{Manifold-adaptive dimension estimation}, at 24th International Conference on Machine Learning 227: 265-272, 2007

    E. Levina and P. J. Bickel, \bold{Maximum Likelihood Estimation of Intrinsic Dimension}, in Advances in Neural Information Processing Systems 17: 777--784, 2005

    D. J. MacKay and Z. Ghahramani, Comments on \dQuote{Maximum Likelihood Estimation of Intrinsic Dimension} by E. Levina and P. Bickel (2004), http://www.inference.phy.cam.ac.uk/mackay/dimension/, 2005

	T. Lin, H. Zha and S. U. Lee, \bold{Riemannian Manifold Learning for Nonlinear Dimensionality Reduction}, at European Conference on Computer Vision, Graz, Austria 9: 44-55, 2006
	
	S-W. Cheng and M-K. Chiu, \bold{Dimension detection via slivers}, at 20th Annual ACM-SIAM Symposium on Discrete Algorithms: 1001-1010, 2009
	
	J. M. Lee and M. Verleysen, \bold{Nonlinear Dimensionality Reduction} (chapter 3), Springer, 2007
}

\examples{
	#generate a swissroll dataset
	n = 300 ; h = 3
	phi = runif(n, min=0, max=2*pi)
	z = runif(n, min=0, max=h)
	sw = cbind( phi*cos(phi), phi*sin(phi), z )

	#estimate dimension
	print(dest_PCA(sw, 20))
	print(dest_unbl(sw, 20))
	print(dest_pett(sw, 20))
	print(dest_fara(sw, 20))
	print(dest_levi(sw, 20))
	print(dest_unbl(sw, 20))
    print(dest_RML(sw, 20, 10))
    print(dest_rgrl(sw, 20, 10))
    print(dest_sliv(sw, 20, 10))
}

