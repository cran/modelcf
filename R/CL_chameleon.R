#update dissimilarities between one cluster (index) and all others
updateDissims = function(data, dissims, clusts, flags, index, alpha) {
    n = nrow(data) ; seqVect = 1:n
    indsRef = seqVect[clusts == clusts[index]] ; lRef = length(indsRef)
    criticSize = 8
    #cluster's edge cut precomputation (see article) :
    EC_ref = 0.0 ; nzRef = 0
    if (lRef >= criticSize) {
        clRef = c()
        if (lRef > 2) clRef = spec_clust("specKM", data[indsRef,], 2, 0, FALSE, getKnn(lRef))
        else clRef = c(1,2)
        locSeq = 1:lRef
        indsI = indsRef[ locSeq[clRef == 1] ]
        indsJ = indsRef[ locSeq[clRef == 2] ]
        srt = sort( as.double(dissims[indsI, indsJ]) )
        card = length(indsI)*length(indsJ)
        nzRef = min(10,card)
        EC_ref = sum(srt[1:nzRef])
    }

    res = rep(Inf, n)
    #now compared to all other clusters :
    M = max(clusts)
    for (i in 1:M) {
        if (i == clusts[index]) next
        inds = seqVect[clusts == i] ; len = length(inds)
        if (len > 0) {
            EC = 0.0 ; nz = 0
            if (len >= criticSize) {
                clu = c()
                if (len > 2) clu = spec_clust("specKM", data[inds,], 2, 0, FALSE, getKnn(len))
                else clu = c(1,2)
                locSeq = 1:len
                indsI = inds[ locSeq[clu == 1] ]
                indsJ = inds[ locSeq[clu == 2] ]
                srt = sort( as.double(dissims[indsI, indsJ]) )
                card = length(indsI)*length(indsJ)
                nz = min(10,card)
                EC = sum(srt[1:nz])
            }

            #fusion of the two clusts ?
            srt = sort( as.double(dissims[indsRef, inds]) )
            card = lRef*len
            normz = min(10,card)
            EC_both = sum(srt[1:normz])
            other_index = inds[ which.max(flags[inds]) ]

            if (nzRef >= criticSize &&  nz >= criticSize) {
                #apply formula :
                RI = (2.0 * EC_both) / (EC_ref + EC)
                RC = ((lRef + len) * EC_both) / (card * (lRef * (EC_ref / nzRef) + len * (EC / nz)))
                res[other_index] = RI * (RC^alpha)
            }
            else {
                #simple update (mean of dissimilarities)
                res[other_index] = EC_both / normz
            }
        }
    }
    return (res)
}

#CHAMELEON clustering : trade-off between closeness and connectivity
chameleon = function(data, dissims, K, alpha) {
    n = nrow(data) ; seqVect = 1:n

    #initialization with a k-means
    clusts = kmeans(data, max(K, min(n-1, floor(n/7))), iter.max=100,nstart=10)$cluster
    flags = rep(FALSE,n) #true if representative element of a cluster
    M = max(clusts)
    for (i in 1:M) {
        indsClust = seqVect[clusts == i]
        flags[indsClust[1]] = TRUE
    }
    diag(dissims) = Inf
    nbDiffClusts = M
    while (nbDiffClusts > K) {

        #fusion one time :
        inds = seqVect[flags]
        fusDissims = dissims[inds,inds]
        toFusion = which.min(fusDissims)
        fusN = length(inds)
        indI = (toFusion-1) %% fusN + 1 #rank
        indJ = floor( (toFusion-1) / fusN ) + 1 #column
        #symmetry => we can switch for the smallest index
        if (indI > indJ) {
            tmp = indI
            indI = indJ
            indJ = tmp
        }

        toFusClust = seqVect[ clusts == clusts[inds[indJ]] ]
        clusts[toFusClust] = clusts[inds[indI]]
        flags[toFusClust] = FALSE
        flags[ seqVect[ clusts == clusts[inds[indI]] ] ] = FALSE
        flags[inds[indI]] = TRUE
        newVectDiss = updateDissims(data, dissims, clusts, flags, inds[indI], alpha)
        dissims[inds[indI],] = newVectDiss
        dissims[,inds[indI]] = newVectDiss
        nbDiffClusts = nbDiffClusts - 1
    }

    return ( reordering(clusts) )
}

