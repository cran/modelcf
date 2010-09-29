## get (dis)similarities from a list of neighborhoods

# symetrize W (already >0 symmetric) "by hand" way :
getSymSims = function(W) {
    n = nrow(W) ; seqVect = 1:n
    for (i in 1:n) {
        #pPos = vector of indices j for which both W[i,j] and W[j,i] are > 0
        pPos = seqVect[W[i,] > 0.0]
        W[i,pPos] = pmin(W[i,pPos],W[pPos,i])
        W[pPos,i] = W[i,pPos]
    }
    return (W)
}

#build similarity matrix W :
gtsimils = function(data, adn, d, knn) {
    n = nrow(data) ; seqVect = 1:n

    #get MUTUAL neighborhoods and their sizes : (=> W is ">0 symmetric")
    NI = getNI(data, adn, d, knn, TRUE)
    knns = sapply(NI, length)

    #get similarity matrix :
    W = matrix(0.0,nrow=n,ncol=n)
    for (i in 1:n) {

        if (knns[i] > 1) {
            localDists = sqrt( colSums( ( t(data[NI[[i]],]) - data[i,] )^2 ) )
            #local sigma :
            firstNN = which.min(localDists >= EPS() )
            sigma2 = (localDists[ knns[i] ]^2 - localDists[firstNN]^2) / log(localDists[ knns[i] ]^2 / localDists[firstNN]^2)
            #recover exceptional errors (shouldn't occur) :
            if (is.nan(sigma2)) sigma2 = 1.0
            #finally compute row of similarity matrix :
            W[ i, NI[[i]] ] = exp( - localDists^2 / sigma2 )
        }

        else if (knns[i] == 1) {
            # only one neighbor
            W[ i, NI[[i]] ] = 1.0
        }
    }

    return ( list("W"=W,"NI"=NI) )
}

#return some technical parameters like stationary distribution, for each graph component
params_rw = function(data, adn, d, knn, symm) {

    #retrieve simils and find number of connected components (each is STRONGLY connected) :
    simils = gtsimils(data, adn, d, knn)
    W = simils$W ; NI = simils$NI
    cc = gt_cxcomps(NI, TRUE) ; nbC = max(cc)
    params = as.list(rep(0,nbC))
    seqVect = 1:nrow(data)

    for (i in 1:nbC) {
        thisComp = seqVect[cc == i]
        nc = length(thisComp)
        if (nc==1) {
            #isolated element
            params[[i]] = list("W"=1.0,"P"=1.0,"ps"=1.0,"sing"=FALSE)
            next
        }

        locW = W[thisComp,thisComp]
        degs = rowSums(locW)
        P = (degs^(-1)) * locW #matrix of transition probabilities
        ps = rep(0.0, nc) #stationary distribution

        #detect periodicity (which would almost never occur in "real world" cases);
        # if periodic, evaluate ps heuristically..
        singular = FALSE
        testP = abs( eigen(P)$values )
        if ( testP[1] - testP[2] < EPS() ) {
            #rough evaluate of ps [which in fact doesn't exist mathematically..] :
            #kind of Cesaro mean ; see Aldous book
            curProb = rep(1.0/nc,nc)
            for (j in 1:nc) curProb = curProb %*% P
            #now average over next nc matrix-vector products :
            for (j in 1:nc) {
                curProb = curProb %*% P
                ps = ps + curProb
            }
            ps = ps / nc
            singular = TRUE
        }
        else {
            #ps = stationary distribution : t(P) ps = ps, so (t(P) - I) ps = 0
            eig = eigen( t(P) - diag(1,nc) )
            ps = abs( Re( eig$vectors[,nc] ) )
            ps = ps / sum(ps)
        }
        ps = as.double(ps)

        if (symm) {
            if (min(ps) >= EPS() ) {
                #symmetrize W accordingly with directed laplacian graph theory :
                tmp = ps * P
                locW = 0.5 * (tmp + t(tmp))
            }
            else {
                #otherwise back to "by hand symmetry" locally :
                locW = getSymSims(locW)
            }
        }
        params[[i]] = list("W"=locW,"P"=P,"ps"=ps,"sing"=singular)
    }

    return ( list("cc"=cc, "params"=params) )
}

#hitting times or commute distances
hitorct = function(data, adn, d, knn, ct, symm, weight) {
    n = nrow(data) ; seqVect = 1:n
    HC = matrix(-Inf,nrow=n,ncol=n)
    prw = params_rw(data,adn,d,knn,symm)
    cc = prw$cc ; nbC = max(cc)

    for (i in 1:nbC) {
        thisComp = seqVect[cc == i]
        nc = length(thisComp)

        if (nc > 1) { #otherwise nothing to do..
            P = prw$params[[i]]$P
            ps = prw$params[[i]]$ps

            Z = matrix() #Fundamental matrix
            if ( ! (prw$params[[i]]$sing) ) Z = mppsinv( diag(1,nc) - t( t(P) + ps ) )
            else {
                #singular case ; we have to use some Cesaro computation..
                Z = P
                for (j in 1:(2*nc)) Z = Z %*% P + P
                Z = Z / (2*nc)
            }

            for (j in 1:nc) {
                if (ps[j] >= EPS() ) HC[ thisComp[-j],thisComp[j] ] = abs( Z[-j,j] - Z[j,j] ) / ps[j]
                else HC[ thisComp[-j],thisComp[j] ] = Inf
            }
            locHC = HC[thisComp,thisComp]
            if (weight) {
                if (ct) HC[thisComp,thisComp] = t( ps * t(locHC) ) + ps * locHC
                else HC[thisComp,thisComp] = t( ps * t(locHC) )
                #..for very singular cases :
                HC[thisComp,thisComp] [ is.nan(HC[thisComp,thisComp]) ] = Inf
            }
            else if (ct) HC[thisComp,thisComp] = locHC + t(locHC)
        }
    }

    M = max( HC[HC < Inf] )
    HC[abs(HC)==Inf] = M + 1.0
    HC = pmax(HC,0.0)
    diag(HC) = 0.0
    return (HC)
}

#get distance matrix from data and similarity : Commute Time
ctdists = function(data, adn, d, knn, symm, weight, sigmo) {
    if (!symm || weight)
        # with hitting times....
        return ( hitorct(data, adn, d, knn, FALSE, symm, weight) )

    #symmetric unweighted case :
    prw = params_rw(data,adn,d,knn,TRUE)
    cc = prw$cc ; nbC = max(cc)
    n = nrow(data) ; seqVect = 1:n
    dists = matrix(-1.0,nrow=n,ncol=n)

    for (i in 1:nbC) {
        thisComp = seqVect[cc == i]
        nc = length(thisComp)

        if (nc > 1) { #otherwise nothing to do..
            #get laplacian :
            W = prw$params[[i]]$W
            L = -W
            diag(L) = rowSums(W)
            invLap = mppsinv(L)

            #if sigmoid-CT-kernel : apply transformation ; we take a==1 (heuristic ?!)
            if (sigmo) {
                sigma = sqrt(var(as.double(invLap)))
                invLap = 1.0 / (1.0 + exp(-invLap / sigma))
            }
            #..and distances
            diagInv = diag(invLap)
            for (i in 1:nc) dists[ thisComp[i],thisComp ] = rep(diagInv[i],nc) + diagInv - 2 * invLap[i,]
        }
    }

    #complete with large distances if not connected
    dists[dists < 0] = max(dists) + 1.0
    diag(dists) = 0.0
    return (dists)
}

