############################################################
## functions to find neighborhoods (normal or LTSA-adaptive)

#get the mutual-knn graph from a list of neighborhoods
#NOTE : some isolated neighborhoods could be empty..
getMutual = function(NI, addl = FALSE) {
    n = length(NI)
    hasI = matrix(FALSE,nrow=n,ncol=n) # [i,j] == i has j as neighb.
    for (i in 1:n) hasI[i,NI[[i]]] = TRUE
    
    if (!addl) hasl = hasI & t(hasI) #deletion of one-sided edges
    else hasl = hasI | t(hasI) #addition of one-sided edges
   
    seqVect = 1:n
    mNI = as.list(rep(0,n))
    for (i in 1:n) mNI[[i]] = seqVect[hasl[i,]]
    return (mNI)
}

#simple k nearest-neighbors graph :
simpleNeighbs = function(data, k, mutual=FALSE) {
    dists = as.matrix(dist(data)) ; diag(dists) = Inf
    n = nrow(data)
    NI = as.list(rep(0,n))
    for (i in 1:n) {
        srt = sort(dists[i,],index.return=TRUE)
        NI[[i]] = srt$ix[1:k]
    }
    if (mutual) return (getMutual(NI))
    return (NI)
}

#aux. function for adaptive kNN :
adaptAux = function(data, d, threshP) {
    n = nrow(data)

    #stage 1: evaluate "relative linearity" at each point
    k = getKnn(n)
    sneighbs = simpleNeighbs(data,k)
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
    k = d
    while (k < n) {
        s = svd(data[ srt$ix[1:(k+1)], ], nu=0, nv=0)
        percent = sum(s$d[1:d]) / sum(s$d)
        if (percent < threshP) break
        k = k + 1
    }

    return ( list("kmax"=k,"percents"=percents) )
}

#a not very good but not too computationaly demanding adaptive neighborhoods :
#WARNING : need an estimation of d
adaptiveNeighbs_basic = function(data,d,kmax=0,mutual=FALSE,threshP = 0.95) {
    n = nrow(data)
    adx = adaptAux(data, d, threshP)
    if (kmax==0) kmax = adx$kmax
    percents = adx$percents

    #now prepare result :
    NI = as.list(rep(0,n))
    sneighbs = simpleNeighbs(data,kmax)
    for (i in 1:n) {
        k = max(1, round(percents[i] * kmax))
        NI[[i]] = sneighbs[[i]] [1:k]
    }
    if (mutual) return (getMutual(NI))
    return (NI)
}

#adaptive neighborhoods, by Zhan et al. (2009)
adaptiveNeighbs1 = function(data,d,kmax=0,mutual=FALSE,threshP = 0.95) {
    n = nrow(data)
    if (kmax==0) kmax = adaptAux(data, d, threshP)$kmax
    kmax = max(d+1, kmax)
	
    #now prepare result :
    NI = as.list(rep(0,n))
    sneighbs = simpleNeighbs(data,kmax)
    for (i in 1:n) {
		NI[[i]] = sneighbs[[i]] [1:d]
		for (j in (d+1):kmax) {
			M = rbind(data[ NI[[i]], ], data[ sneighbs[[i]] [j], ])
			s = svd(M,nu=0,nv=0)
			if (sum(s$d[1:d]) / sum(s$d) < threshP) break
			NI[[i]] = c(NI[[i]], sneighbs[[i]] [j])
		}
    }
    if (mutual) return (getMutual(NI))
    return (NI)
}

#adaptive neighborhoods (from "Adaptive manifold Learning", adaptive LTSA algorithm)
#WARNING : need an estimation of d
adaptiveNeighbs2 = function(data,d,kmax=0,mutual=FALSE, eta = 0.05, expand=FALSE) {
    n = nrow(data)
    if (kmax==0) kmax = max(d+1, getKnn(n))
    kmin = min( 3, kmax-d, max(1, n %/% 30) )
    ratio = matrix(0,nrow=n,ncol=kmax)
    NI_full = simpleNeighbs(data,kmax)
    for (i in 1:n) NI_full[[i]] = c(i, NI_full[[i]])
    NI = as.list(rep(0,n))
    knns = rep(0,n)
    k = min(n - 1, kmax)
    eps_norms = as.list(rep(0,n))
    coefs_norms = as.list(rep(0,n))

    #NEIGHBORHOOD CONTRACTION
    isOver = rep(FALSE,n)
    oldIsOver = rep(TRUE, n)
    while ( sum(isOver==oldIsOver) < n) {
        oldIsOver = isOver
        for (i in 1:n) {

            if (!isOver[i]) {
                mean = colMeans(data[ NI_full[[i]] [1:(k+1)], ])
                ctrData = t( t(data[ NI_full[[i]] [1:(k+1)], ]) - mean )
                s = svd(ctrData,nu=0,nv=d) ; V = as.matrix(s$v)
                coefs = ctrData %*% V
                #test eta-dependant : if reconstructions good enough, stop
                eps_norm = sqrt( sum( ( ctrData - coefs %*% t(V) )^2 ) )
                coefs_norm = sqrt( sum(coefs^2) )
                if ( eps_norm <= eta * coefs_norm ) {
                    #terminate with current value of knn for i :
                    NI[[i]] = NI_full[[i]] [2:(k+1)]
                    knns[i] = k
                    isOver[i] = TRUE
                }
                else ratio[i,k] = eps_norm / coefs_norm
            }
        }
        if (k <= d + kmin) break
        k = k - 1
    }

    seqVect = 1:n
    #fill the remaining neighborhoods if needed :
    for (i in seqVect[ !isOver ]) {
        knns[i] = d+kmin-1 + which.min( ratio[i,(d+kmin):kmax] )
        NI[[i]] = NI_full[[i]] [2:(knns[i]+1)]
    }

    if (expand) {
        #NEIGHBORHOOD EXPANSION (costly ! and seems not so useful)
        for (i in 1:n) {
            if (knns[i] == kmax) next
            toAdd = c()
            for (j in (knns[i]+1):kmax) {
                mean = colMeans(data[ NI_full[[i]] [1:(j+1)], ])
                ctrData = t( t(data[ NI_full[[i]] [1:(j+1)], ]) - mean )
                s = svd(ctrData,nu=0,nv=d) ; V = as.matrix(s$v)
                coef_j = ctrData %*% V
                eps_norm_j = sqrt( sum( (ctrData - coef_j %*% t(V) )^2 ) )
                if ( eps_norm_j <= sqrt( sum(coef_j^2) ) ) toAdd = c(toAdd,j)
            }
            NI[[i]] = c(NI[[i]], NI_full[[i]] [toAdd])
        }
    }

    if (mutual) return (getMutual(NI))
    return (NI)
}

#find neighborhoods, simple or adaptive way
getNI = function(data, adn, d, k, mutual=FALSE, threshP = 0.95, eta = 0.05, expand=FALSE) {
    NI = list()
    if (k==0) k = getKnn( nrow(data) )
    if (adn != "none") {
        n = nrow(data)
        adx = adaptAux(data, d, threshP)
        kmax = adx$kmax
        NI = list()
        if (adn=="adbas") NI = adaptiveNeighbs_basic(data, d, kmax, mutual, threshP)
        else if (adn=="ad1") NI = adaptiveNeighbs1(data, d, kmax, mutual, threshP)
        else if (adn=="ad2") NI = adaptiveNeighbs2(data, d, kmax, mutual, eta, expand)
    }
    else NI = simpleNeighbs(data, k, mutual)
    return (NI)
}

#neighborhoods "RML style" :
neighbs_RML = function(data, rgdists, kmin, kmax, tsoft=0.1) {
    n = nrow(data)
    prodTab = n * kmax
    tabK = .C("estimK_RML",data=as.double(t(data)),rgdists=as.double(t(rgdists)),n=as.integer(n),
				m=as.integer(ncol(data)),kmin=as.integer(kmin),kmax=as.integer(kmax),tsoft=as.double(tsoft),
                tabK=integer(prodTab), PACKAGE=pkgnm())$tabK

    NI = as.list(rep(0,n))
    for (i in 1:n) {
        NI[[i]] = tabK[ ((i-1)*kmax+1) : (i*kmax) ]
        NI[[i]] = NI[[i]] [NI[[i]] > 0]
    }
    return (NI)
}

