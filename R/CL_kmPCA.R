#function by Chiou and Li : PCA, generalizing classical k-means for functional case
#can also be applied to reduced coordinates (with a [very] small loss)
#will be very (very) slow if simplif==FALSE on high-dim. data (original algorithm)
#simplif==TRUE will cut the leave-one-out estimate of PCA basis
km_PCA = function(data, K, d, simplif=TRUE, maxiter=50) {
    # initializations I
    n = nrow(data) ; m = ncol(data) ; seqVect = 1:n
    clusts = as.list(rep(0,K))
    locBasis = as.list(rep(0,K))
    ctrs = matrix(nrow=K, ncol=m)

    #INITIALIZATIONS II
    km = kmeans(data,K,iter.max=100,nstart=10)
    ctrs = km$centers
    #estimate local bases :
    for (i in 1:K) {
        clusts[[i]] = seqVect[km$cluster == i]
        ctrData = t( t(data[clusts[[i]],]) - colMeans(data[clusts[[i]],]) )
        s = svd( ctrData, nu=0, nv = min(d, length(clusts[[i]]) ) )
        locBasis[[i]] = as.matrix(s$v)
    }

    ## main loop :
    errs = matrix(nrow=n, ncol=K)
    lastCl = rep(0,n) ; newCl = rep(0,n)
    counter = 1
    while (counter <= maxiter) {

        errs[,] = Inf
        for (i in 1:K) {
            #compute errors for every point [with special case if points are within the current cluster] :
            normalComp = c()
            if (simplif) normalComp = seqVect
            else normalComp = seqVect[ - clusts[[i]] ]

            ctrData = t( t(data[normalComp,]) - ctrs[i,] )
            errFunc = ctrData - (ctrData %*% locBasis[[i]]) %*% t(locBasis[[i]])
            errs[normalComp,i] = sqrt( rowSums(errFunc^2) / m )

            if (!simplif) {
                for (j in clusts[[i]]) {
                    inds = clusts[[i]] [-j]
                    if (length(inds) == 0) break

                    #compute "leave-one-out" mean and basis :
                    loo_mean = colMeans(data[inds,])
                    ctrData = t( t(data[inds,]) - loo_mean )
                    s = svd( ctrData, nu=0, nv=min(d,length(inds)) )
                    loo_basis = as.matrix(s$v)

                    #compute best approximation error :
                    delta = data[j,] - loo_mean
                    errFunc = delta - (delta %*% loo_basis) %*% t(loo_basis)
                    errs[j,i] = sqrt(sum(errFunc^2) / m)
                }
            }
        }

        #estimate new clusters memberships
        for (i in 1:K) {
            lastCl[ clusts[[i]] ] = i #take a "print" of last clusters
            clusts[[i]] = 0
        }
        for (i in 1:n) {
            ind = which.min(errs[i,])
            if (clusts[[ind]] [1] == 0) clusts[[ind]] = i
            else clusts[[ind]] = c(clusts[[ind]], i)
        }

        #..and new bases
        infK = FALSE
        for (i in 1:K) {
            if (clusts[[i]] [1] == 0) {
                infK = TRUE
                break
            }
            ctrs[i,] = colMeans(data[clusts[[i]],])
            ctrData = t( t(data[ clusts[[i]], ]) - ctrs[i,] )
            s = svd( ctrData, nu=0, nv=min(d,length(clusts[[i]])) )
            locBasis[[i]] = as.matrix(s$v)
        }

        #emergency case : fewer clusters than K, return simple k-means :
        if (infK) return (km$cluster)

        #compare new partition to old one (simple matching : ordering is preserved) :
        for (i in 1:K) newCl[ clusts[[i]] ] = i
        ratio = sum( newCl - lastCl == 0) / n #rate of "well classified"
        if (ratio >= 0.95) return (newCl)

        counter = counter+1
    } ## end main loop

    return (newCl)
}

