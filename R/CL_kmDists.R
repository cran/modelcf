# k-means based on a distance matrix
km_dists = function(distm, K, nstart=10, maxiter=100) {
    # initializations :
    n = nrow(distm)
    clusts = as.list(rep(0,K))
    bestClusts = 0 ; bestDistor = 0.0

    for (repet in 1:nstart) {
        ctrs = sort(sample(1:n,K)) # centers initialization
        oldCtrs = rep(0,K)
        counter = 1

        ## main loop :
        # while old and "new" centers differ..
        while (sum(ctrs!=oldCtrs) > 0 && counter<=maxiter) {

            # clusters reinitialization :
            for (i in 1:K) clusts[[i]] = 0

            #estimate clusters belongings
            for (i in 1:n) {
                indMin = which.min(distm[i,ctrs])
                if (clusts[[indMin]] [1] == 0) clusts[[indMin]] = i
                else clusts[[indMin]] = c(clusts[[indMin]],i)
            }

            oldCtrs = ctrs
            # recompute centers :
            for (i in 1:K) {
                L = length(clusts[[i]])
                if (L > 1) ctrs[i] = clusts[[i]] [ which.min(rowSums(distm[ clusts[[i]],clusts[[i]] ])) ]
                else ctrs[i] = clusts[[i]] [1] #L >=1 because contains at least the central point
            }
            ctrs=sort(ctrs)
            counter = counter+1
        } ## end main loop

        # ..finally compute distorsions :
        evals = rep(0.0,K)
        for (i in 1:n) {
            indMin = which.min(distm[i,ctrs])
            evals[indMin] = evals[indMin] + min(distm[i,ctrs])
        }
        distor = sum(evals)
        if (bestDistor==0.0 || distor < bestDistor) {
            bestClusts = clusts
            bestDistor = distor
        }
    }

    # and return clusters
    res = rep(0,n)
    for (i in 1:K) {
        for (j in 1:length(bestClusts[[i]]))
            res[ bestClusts[[i]] [j] ] = i
    }
    return (res)
}

