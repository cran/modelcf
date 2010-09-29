##SPECTRAL CLUSTERING FOR (possibly) DISCONNECTED GRAPHS :

#"custom" spectral embedding
spec_emb = function(data, K, d, adn, knn) {
    n = nrow(data) ; seqVect = 1:n
    embedding = matrix(0.0, nrow=n, ncol=K);
    #get parameters for each component (W, P, ps)
    prw = params_rw(data,adn,d,knn,TRUE)
    cc = prw$cc ; nbC = max(cc)

    for (i in 1:nbC) {
        thisComp = seqVect[cc == i]
        nc = length(thisComp)
        if (nc==1) {
            #isolated element : do nothing
            next
        }

        P = prw$params[[i]]$P
        ps = prw$params[[i]]$ps
        locW = prw$params[[i]]$W

        # and get local laplacian :
        D = rowSums(locW)
        L = -locW ; diag(L) = D
        D[D < EPS()] = Inf ; D = 1.0/D
        D = sqrt(D)
        L = t( D * t(D * L) )
        eig = eigen(L,symmetric = TRUE)

        #finally get (local) embedding from the eigenvectors
        embedding[thisComp,] = eig$vectors[, nc - ( 1 : min(K, nc-1) ) ]
    }

    return ( list("embed"=embedding, "cc"=cc) );
}

#..then spectral clustering
spec_clust = function(method, data, K, d, adn, knn) {
    embcc = spec_emb(data, K, d, adn, knn)
    embed = embcc$embed ; cc = embcc$cc
    seqVect = 1:nrow(data)
    nbC = max(cc);

    # first check if more components than clusters asked :
    if (nbC >= K) {
        # "easy", juste merge some clusters :
        return ( mergeToK(data, cc, K) );
    }

    #otherwise, we have to cluster somewhere..
    #we choose the larger one (and repeat if needed !)
    curNbC = nbC
    while (curNbC < K) {
        #evaluate size of larger clust VS next larger one :
        maxInd1 = 0 ; maxInd2 = 0
        maxSize1 = 0 ; maxSize2 = 0
        for (i in 1:curNbC) {
            curSize = sum(cc==i)
            if (curSize > maxSize1) {
                maxSize2 = maxSize1
                maxSize1 = curSize;
                maxInd2 = maxInd1
                maxInd1 = i
            }
            else if (curSize > maxSize2) {
                maxSize2 = curSize
                maxInd2 = i
            }
        }

        #cut the larger one in max(2, min(floor(maxSize1/maxSize2), K-curNbC+1) ) parts :
        cutIn = max(2, min(floor(maxSize1/maxSize2), K-curNbC+1) )
        inds = seqVect[cc==maxInd1]
        clusts = c()
        if (method=="specH") clusts = phclust(dist(embed[inds,1:cutIn]), cutIn)
        else if (method=="specKM") clusts = (kmeans(embed[inds,1:cutIn], cutIn, iter.max=100, nstart=10))$cluster
        #split cc accordingly :
        cc[inds] = max(cc) + clusts - 1
        curNbC = curNbC + max(clusts) - 1;
    }

    return ( reordering(cc) )
}

