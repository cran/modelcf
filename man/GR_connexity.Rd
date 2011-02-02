\name{connexity}
\alias{gt_cxcomps}
\alias{testConnexity}
\title{Functions around graph conenxity}

\description{
    Internal use only ; should not be called by the final user.

    \code{gt_cxcomps} gets the connected components based on neighborhoods.

    \code{testConnexity} returns neighborhoods that assure graph (weak) connexity. If \code{NI} does not lead to connexity, kNN-graph is built instead.
}

\usage{
    gt_cxcomps(NI, ctype=FALSE, k=0)

    testConnexity(data, NI, k)
}

\arguments{
    \item{data}{matrix of n vectors in rows ; \code{data[i,]} is the i-th m-dimensional vector}

    \item{NI}{list of (graph) neighborhoods}

    \item{k}{fixed number of neighbors at each point ; used only if \code{adn == FALSE}}

    \item{ctype}{(weak) connexity type, TRUE for mutual-kNN graph (clustering case), FALSE for general graphs (dimensionality reduction)}
}

\value{
    \code{gt_cxcomps} returns a vector like (1 1 1 2 2 2 1 1 3 3 3) describing connected components.

    \code{testConnexity} returns the smallest list of neighborhoods that ensure connexity.
}

