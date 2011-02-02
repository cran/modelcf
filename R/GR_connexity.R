#return the number of connected components of the kNN-graph :
#only "weak connexity" test because we use mutual neighborhoods for clustering,
#and don't care the edges directions for dimensionality reduction
#ctype==TRUE for weak def., mutual-kNN graph (==undirected, clustering),
#ctype==FALSE for weak def., any graph (==undirected, red. dim.)
gt_cxcomps = function(NI, ctype=FALSE, k=0) {
    n = length(NI)
    binarySim = matrix(0,nrow=n,ncol=n)
    if (k > 0) { #for getLimitConnex :
        for (i in 1:n) binarySim[ NI[[i]][1:k], i] = 1
    }
    else { #all other cases (adaptive or not)
        for (i in 1:n) binarySim[ NI[[i]], i] = 1
    }

    cc = .C("getConnex",binMat=as.double(binarySim),n=as.integer(n),
            ctype=as.integer(ctype), cc=integer(n), PACKAGE=pkgnm())$cc
    return ( reordering(cc) )
}

#get "limit" connexity, when smaller k lead to fragmented graph :
#WARNING : n >= 4..
getLimitConnex = function(data) {
    n=nrow(data)
    left = 1; right = n-1
    midK=0; oldMidK = 0

    # compute the local neighborhoods "at every scale" :
    tabK_ref = simpleNeighbs(data,n-1)

    #dichotomic search :
    repeat {
        midK = ceiling( (left+right)/2 )
        cc = max(gt_cxcomps(tabK_ref, FALSE, midK))
        if (cc > 1) left = midK
        else right = midK

        if (midK == oldMidK) break
        oldMidK = midK
    }
    return (midK)
}

#return neighborhoods that assure graph (weak) connexity
testConnexity = function(data, NI, k) {
    cc = max(gt_cxcomps(NI))
    if (cc > 1) {
        #if not connex, back to simple neighborhoods (for dim. red.)
        limit_kNN = getLimitConnex(data)
        k = max(k, limit_kNN)
        NI = simpleNeighbs(data,k)
    }
    return (NI)
}

