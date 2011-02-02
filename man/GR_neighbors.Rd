\name{neighbors}
\alias{getNI}
\alias{neighbs_RML}
\title{Top-level neighborhoods functions}

\description{
    Internal use only ; should not be called by the final user.

    \code{getNI} gets neighborhoods (designed for any algorithm), adaptively or not.

    \code{neighbs_RML} gets neighborhoods designed for Riemannian Manifold Learning algorithm (see \code{\link{RML}}).
}

\usage{
    getNI(data, adn, d, k, mutual=FALSE, threshP = 0.95, eta = 0.05,
    expand=FALSE)

    neighbs_RML(data, rgdists, kmin, kmax, tsoft=0.1)
}

\arguments{
    \item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}
    
    \item{rgdists}{rough approximations of graph distances (only for \code{RML})}

    \item{kmin}{minimum number of neighbors at each point}

    \item{kmax}{maximum number of neighbors at each point}

    \item{adn}{string for adapted point-varying neighborhoods. "none" for no adaptivity, "adbas" for simple local PCA based neighborhoods (see code), "ad1" for the Zhan et al. method, and "ad2" for Wang et al. method. In short, the more linear data is around \code{x}, the more \code{x} has neighbors}

    \item{d}{estimated data dimensionality ; useful only if \code{adn} is set}

    \item{k}{fixed number of neighbors at each point ; used only if \code{adn == FALSE}}

    \item{tsoft}{tolerance factor for the visibility graph computation (between 0 and 1 ; should be close to 0)}

    \item{mutual}{boolean for computing mutual graph neighborhoods}

    \item{threshP}{threshold percentage to estimate number of neighbors for a locally linear surface}

    \item{eta}{threshold percentage used for neighborhood contraction}

    \item{expand}{boolean for enable the final expansion step in the Wang et al. algorithm}
}

\details{
    \code{neighbs_RML} computes the visibility graph as described in the article of Lin et al. See this paper for further explanations.
}

\value{
    A list of neighborhoods, describing a graph.
}

\references{
    T. Lin, H. Zha and S. U. Lee, \bold{Riemannian Manifold Learning for Nonlinear Dimensionality Reduction}, at European Conference on Computer Vision, Graz, Austria 9: 44-55, 2006

    J. Wang, Z. Zhang and H. Zha, \bold{Adaptive Manifold Learning}, in Advances in Neural Information Processing Systems 17: 1473-1480, 2005

    Y. Zhan, J. Yin, X. Liu and G. Zhang, \bold{Adaptive Neighborhood Select Based on Local Linearity for Nonlinear Dimensionality Reduction}, at International Symposium on Advances in Computation and Intelligence, Huangshi, China 5821: 337-348, 2009
}

