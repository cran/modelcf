############################################################
## functions to find neighborhoods (normal or LTSA-adaptive)

#get the mutual-knn graph from a list of neighborhoods
#NOTE : some isolated neighborhoods could be empty..
getMutual = function(NI) {
    n = length(NI)
    hasI = matrix(FALSE,nrow=n,ncol=n) # [i,j] == i has j as neighb.
    for (i in 1:n) hasI[i,NI[[i]]] = TRUE
    hasI = hasI & t(hasI) #deletion of one-sided edges
    seqVect = 1:n
    mNI = list()
    for (i in 1:n) mNI[[i]] = seqVect[hasI[i,]]
    return (mNI)
}

#simple k nearest-neighbors graph :
simpleNeighbs = function(data, knn, mutual=FALSE) {
    dists = as.matrix(dist(data)) ; diag(dists) = Inf
    n = nrow(data)
    NI = as.list(rep(0,n))
    for (i in 1:n) {
        srt = sort(dists[i,],index.return=TRUE)
        NI[[i]] = srt$ix[1:knn]
    }
    if (mutual) return (getMutual(NI))
    return (NI)
}

#aux. function for adaptive kNN :
adaptAux = function(data, d, threshP) {
    n = nrow(data)

    #stage 1: evaluate "relative linearity" at each point
    kNN = min(n-1, max(10,3*d))
    sneighbs = simpleNeighbs(data,kNN)
    percents = rep(0.0,n)
    for (i in 1:n) {
        sneighbs[[i]] = c(i,sneighbs[[i]])
        s = svd(data[ sneighbs[[i]], ], nu=0, nv=0)
        percents[i] = sum(s$d[1:d]) / sum(s$d)
    }
    percents = percents / max(percents)
    testInd = which.max(percents)

    #stage 2: determine kNN at the "best linear" point, then infer all others kNN
    dists = sqrt( colSums( (t(data) - data[testInd,])^2 ) )
    srt = sort(dists,index.return=TRUE)
    percent = 1.0
    kNN = d
    while (kNN < n) {
        s = svd(data[ srt$ix[1:(kNN+1)], ], nu=0, nv=0)
        percent = sum(s$d[1:d]) / sum(s$d)
        if (percent < threshP) break
        kNN = kNN + 1
    }

    return ( list("knnmax"=kNN,"percents"=percents) )
}

#a not very good but not too computationaly demanding adaptive neighborhoods :
#WARNING : need an estimation of d
adaptiveNeighbs1 = function(data,d,mutual=FALSE,threshP = 0.95) {
    n = nrow(data)
    adx = adaptAux(data, d, threshP)
    knnmax = adx$knnmax
    percents = adx$percents

    #now prepare result :
    NI = as.list(rep(0,n))
    sneighbs = simpleNeighbs(data,knnmax)
    for (i in 1:n) {
        knn = max(1, round(percents[i] * knnmax))
        NI[[i]] = sneighbs[[i]] [1:knn]
    }
    if (mutual) return (getMutual(NI))
    return (NI)
}

#adaptive neighborhoods (from "Adaptive manifold Learning", adaptive LTSA algorithm)
#WARNING : need an estimation of d
adaptiveNeighbs2 = function(data,d,knnmax,mutual=FALSE, threshP=0.95, eta = 0.05, expand=FALSE) {
    n = nrow(data)
    knnmax = max(d+1, knnmax)
    knnmin = min( 3, knnmax-d, max(1, n %/% 30) )
    ratio = matrix(0,nrow=n,ncol=knnmax)
    NI_full = simpleNeighbs(data,knnmax)
    for (i in 1:n) NI_full[[i]] = c(i, NI_full[[i]])
    NI = as.list(rep(0,n))
    knns = rep(0,n)
    knn = min(n - 1, knnmax)
    eps_norms = as.list(rep(0,n))
    coefs_norms = as.list(rep(0,n))

    #NEIGHBORHOOD CONTRACTION
    isOver = rep(FALSE,n)
    oldIsOver = rep(TRUE, n)
    while ( sum(isOver==oldIsOver) < n) {
        oldIsOver = isOver
        for (i in 1:n) {

            if (!isOver[i]) {
                mean = colMeans(data[ NI_full[[i]] [1:(knn+1)], ])
                ctrData = t( t(data[ NI_full[[i]] [1:(knn+1)], ]) - mean )
                s = svd(ctrData,nu=0,nv=d) ; V = as.matrix(s$v)
                coefs = ctrData %*% V
                #test eta-dependant : if reconstructions good enough, stop
                eps_norm = sqrt( sum( ( ctrData - coefs %*% t(V) )^2 ) )
                coefs_norm = sqrt( sum(coefs^2) )
                if ( eps_norm <= eta * coefs_norm ) {
                    #terminate with current value of knn for i :
                    NI[[i]] = NI_full[[i]] [2:(knn+1)]
                    knns[i] = knn
                    isOver[i] = TRUE
                }
                else ratio[i,knn] = eps_norm / coefs_norm
            }
        }
        if (knn <= d + knnmin) break
        knn = knn - 1
    }

    seqVect = 1:n
    #fill the remaining neighborhoods if needed :
    for (i in seqVect[ !isOver ]) {
        knns[i] = d+knnmin-1 + which.min( ratio[i,(d+knnmin):knnmax] )
        NI[[i]] = NI_full[[i]] [2:(knns[i]+1)]
    }

    if (expand) {
        #NEIGHBORHOOD EXPANSION (costly ! and seems not so useful)
        for (i in 1:n) {
            if (knns[i] == knnmax) next
            toAdd = c()
            for (j in (knns[i]+1):knnmax) {
                mean = colMeans(data[ NI_full[[i]] [1:(j+1)], ])
                ctrData = t( t(data[ NI_full[[i]] [1:(j+1)], ]) - mean )
                s = svd(ctrData,nu=0,nv=d) ; V = as.matrix(s$v)
                coef_j = data[NI_full[[i]] [j+1],] %*% V
                eps_norm_j = sqrt( sum( (data[NI_full[[i]] [j+1],] - coef_j %*% t(V) )^2 ) )
                if ( eps_norm_j <= sqrt( sum(coef_j^2) ) ) toAdd = c(toAdd,j)
            }
            NI[[i]] = c(NI[[i]], NI_full[[i]] [toAdd])
        }
    }

    if (mutual) return (getMutual(NI))
    return (NI)
}

#find neighborhoods, simple or adaptive way
getNI = function(data, adn, d, knn, mutual=FALSE, threshP = 0.95, eta = 0.05) {
    NI = list()
    if (knn==0) knn = getKnn( nrow(data) )
    if (adn) {
        n = nrow(data)
        adx = adaptAux(data, d, threshP)
        knnmax = adx$knnmax
        NI = adaptiveNeighbs2(data, d, knnmax, mutual, threshP, eta)
    }
    else NI = simpleNeighbs(data, knn, mutual)
    return (NI)
}

#neighborhoods "RML style" :
neighbs_RML = function(data, knnmin, knnmax, tsoft=0.1) {
    n = nrow(data)
    prodTab = n * knnmax
    tabK = .C("estimK_RML",data=as.double(t(data)),n=as.integer(n),m=as.integer(ncol(data)),
                knnmin=as.integer(knnmin),knnmax=as.integer(knnmax),tsoft=as.double(tsoft),
                tabK=integer(prodTab), PACKAGE=pkgnm())$tabK

    NI = as.list(rep(0,n))
    for (i in 1:n) {
        NI[[i]] = tabK[ ((i-1)*knnmax+1) : (i*knnmax) ]
        NI[[i]] = NI[[i]] [NI[[i]] > 0]
    }
    return (NI)
}

