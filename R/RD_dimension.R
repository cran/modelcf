########################################################
## some methods for dimension estimation :

#local PCA method using a (nonoptimal) covering of the data
dest_PCA = function(data, k, thvar=0.01) {
    n = nrow(data)
    nruns = max(2, 2*floor(n/k))
    NI = simpleNeighbs(data, k)

    #estimate every local dimension, then average
    dim = 0
    picked = rep(FALSE,n)
    refVect = 1:n
    counter = 0
    for (i in 1:nruns) {
        pt = sample(refVect[picked==FALSE],1)
        ptVois = c(pt,NI[[pt]])
        picked[ptVois] = TRUE
        locData = t( t(data[ptVois,]) - colMeans(data[ptVois,]) )
        s = svd(locData,nu=0,nv=0)
        locdim = 1
        locThresh = thvar * s$d[1]
        maxDim = min(nrow(locData),ncol(locData))
        while(locdim <= maxDim && s$d[locdim] >= locThresh)
            locdim = locdim + 1
        dim = dim + locdim - 1
        counter = counter + 1
        if (sum(picked) == n) break
    }
    return ( round(dim/counter) )
}

#local PCA method from Bruske & Sommer, 1998
dest_clust = function(data, nclusts, thvar=0.01) {
    #computation of the "windows" (using neighborhoods) :
    n = nrow(data)
    k = kmeans(data, nclusts, iter.max=100, nstart=10)

    #estimate every local dimension, then average
    dim = 0
    for (i in 1:nclusts) {
        locData = data[k$cluster==i,]
	if (is.null(dim(locData))) {
		dim = dim + 1
		break
	}
        s = svd(locData,nu=0,nv=0)
        locdim = 0
        locThresh = thvar * s$d[1]
        maxDim = min(nrow(locData),ncol(locData))
        while(locdim < maxDim && s$d[locdim+1] >= locThresh)
            locdim = locdim + 1
        dim = dim + locdim
    }
    return ( round(dim/nclusts) )
}

#slope of linear regression Y~X
slope = function(X, Y) {
    n = length(X)
    SX = sum(X)
    SY = sum(Y)
    S2X = sum(X^2)
    SXY = sum(X * Y)
    return ( (SX*SY - n*SXY) / (SX^2 - n*S2X) )
}
#aux function to estimate dimension using linear regression :
linDim = function(data, kmax, d) {
    n = nrow(data)
    G = c()
    rgK = 1:kmax
    if (d > 0) G = (rgK^(1.0/d) * gamma(rgK)) / gamma(rgK + 1.0/d)
    else G = rep(0.0, kmax)
    mdists = as.matrix(dist(data)) ; diag(mdists) = Inf
    for (i in 1:n) mdists[i,] = sort(mdists[i,])
    mdists = colMeans(mdists)
    invD = slope(log(rgK), G+log(mdists[rgK]))
    return ( round(1.0/invD) )
}
#An intrinsic dimensionality estimator from near-neighbor information, Pettis et al. (1979)
dest_pett = function(data, kmax) {
    oldD = 0
    d = linDim(data, kmax, oldD)
    maxIter = 50
    i = 0
    while (i < maxIter && abs(d-oldD) > 0) {
        oldD = d
        d = linDim(data, kmax, oldD)
        i = i + 1
    }
    return (d)
}

#Maximum Likelihood Estimation of Intrinsic Dimension, Levina & Bickel (2005)
dest_levi = function(data, k) {
    dists = as.matrix(dist(data))
    diag(dists) = Inf
    n = nrow(data)
    k = max(2,k)

    #compute all local dimensions and average them
    dim = 0.0
    for (i in 1:n) {
		locDim = 0.0
		srt = sort(as.double(dists[i,]))
		for (j in 1:(k-1))
			locDim = locDim + log(srt[k]/srt[j])
		dim = dim + (k-1)/locDim
	}

	return (round (dim/n) )
}

#debiasing of the last method, by MacKay and Ghahramani (2005) :
dest_unbl = function(data, k) {
    dists = as.matrix(dist(data))
    diag(dists) = Inf
    n = nrow(data)
    k = max(2,k)

    #compute all local inverse dimensions and average them
    invDim = 0.0
    for (i in 1:n) {
		locDim = 0.0
		srt = sort(as.double(dists[i,]))
		for (j in 1:(k-1))
			locDim = locDim + log(srt[k]/srt[j])
		invDim = invDim + locDim/(k-1)
	}

	return (round (n/invDim) )
}

#Manifold-Adaptive Dimension Estimation, Farahmand et al. (2006)
dest_fara = function(data, k) {
    dists = as.matrix(dist(data))
    diag(dists) = Inf
    n = nrow(data)
    k = max(2,k)

    #compute all local dimensions and average them
    dim = 0
    refVect = 1:n
    for (i in 1:n) {
        classI = rank(dists[i,],ties.method="random")
        j1 = refVect[classI == k%/%2]
        j2 = refVect[classI == k]
        locdim = log(2) / (log(dists[i,j2] / dists[i,j1]))
        dim = dim + locdim
    }
    return ( round(dim/n) )
}

#"stochastic" RML method
dest_RML = function(data, kmax, N=10, tsoft=0.0) {
	n = nrow(data)
	simpCount = rep(0,n)
    infini = max(dist(data)) + 1.0

    #initialize distances :
    NI = neighbs_RML(data, 1, kmax, tsoft)
    NI = getMutual(NI)
    dists = matrix(infini,nrow=n,ncol=n)
    for (i in 1:n) {
        if (length(NI[[i]]) > 1)
            dists[ i, NI[[i]] ] = sqrt( colSums( ( t(data[NI[[i]],]) - data[i,] )^2 ) )
        else dists[i, NI[[i]] ] = sqrt( sum( (data[NI[[i]],] - data[i,] )^2 ) )
    }

    for (i in 1:N) {
        #greedy heuristic, very similar to hierarchical clustering
        simpCount = simpCount + .C("rand_hclust", dissims=as.double(t(dists)), n=as.integer(n),
                                 infini=as.double(infini), stats=integer(n), PACKAGE=pkgnm())$stats
    }
    simpCount = simpCount / N

    #regarder simpCount, l'element le + haut significatif = la dimension ; ex 12 8 4 1 0 ==> 2D...
	for (d in 3:n) {
		if (simpCount[d] < 1.0) return (d-2)
	}
	return (n-1)
}

#regularization attempt
dest_rgrl = function(data, kmax, N=10, alpha=3) {
	n = nrow(data)
	simpCount = rep(0,n)
    infini = max(dist(data)) + 1.0

    tsvals = seq(from=0.0, to=0.9, by=0.1)
    ponds = rep(0.0,length(tsvals))
    curind = 1
    for (tsoft in tsvals) {
        #initialize distances :
        NI = neighbs_RML(data, 1, kmax, tsoft)
        NI = getMutual(NI) #avoid noisy simplices
        dists = matrix(infini,nrow=n,ncol=n)
        for (i in 1:n) {
            if (length(NI[[i]]) > 1)
                dists[ i, NI[[i]] ] = sqrt( colSums( ( t(data[NI[[i]],]) - data[i,] )^2 ) )
            else if (length(NI[[i]]) == 1)
                dists[i, NI[[i]] ] = sqrt( sum( (data[NI[[i]],] - data[i,] )^2 ) )
        }

        ponds[curind] = (1.0 - tsoft)^alpha
        for (i in 1:N) {
            #greedy heuristic, very similar to hierarchical clustering
            simpCount = simpCount + ponds[curind] * .C("rand_hclust", dissims=as.double(t(dists)), n=as.integer(n),
                                                     infini=as.double(infini), stats=integer(n), PACKAGE=pkgnm())$stats
        }
        curind = curind + 1
    }
    simpCount = simpCount / (N*sum(ponds))

	for (d in 3:n) {
		if (simpCount[d] < 1.0) return (d-2)
	}
	return (n-1)
}

#generate random data according to a radii distibution
genRdata = function(n, d, densit) {
	#on a sphere :
	coefs = matrix(runif(n*(d+1),min=0.0,max=2*pi),nrow=n,ncol=d+1)
	rData = matrix(nrow=n, ncol=d+1)
	for (i in 1:(d+1)) {
		prod = rep(1.0, n)
		if (i > 1) {
			for (j in 1:(i-1))
				prod = prod * sin(coefs[,j])
		}
		if (i <= d) prod = prod * cos(coefs[,i+1])
		rData[,i] = prod
	}

	#spread according to density :
	radii = c()
	M = max(densit$y)
	L = length(densit$x)
	lvect = 1:L
	while (length(radii) < n) {
		abscisses = sample(lvect, n, replace=TRUE)
		ordonnees = runif(n, min=0.0, max=M)
		accepted = (ordonnees <= densit$y[abscisses])
		if (sum(accepted) > 0) {
			accepted = lvect[ abscisses[accepted] ]
			radii = c(radii, densit$x[accepted])
		}
	}
	return (radii[1:n] * rData)
}
#compute the volume of a simplex from distances matrix
simpVol = function(data) {
	n = nrow(data)
	dists2 = ( as.matrix( dist(data)^2 ) )
	M = cbind( rbind(dists2, rep(1, n)), c(rep(1,n),0) )
	return ( sqrt(abs(det(M))) / (sqrt(2^(n-1)) * factorial(n-1)) )
}
#sample N random volumes accordin to a radii density
sampleVols = function(N, d, densit) {
	res = c()
	for (i in 1:N) {
		rdata = genRdata(d+1,d, densit)		
		res = c(res, simpVol(rdata) )
	}
	return (res)
}
#trial to estimate dimension based on simplices volumes / slivers detection :
dest_sliv = function(data, k, N=10000, thtest=0.05) {
    NI = simpleNeighbs(data, k)
	#distances :
	dists = c()
	n = nrow(data)
	D = ncol(data)
    for (i in 1:n) {
		NI[[i]] = c(i, NI[[i]])
		dists = c( dists, dist( data[NI[[i]],] ) )
	}
	#estimated distances density
	densit = density(dists, kernel = "epanechnikov", from=min(dists), to=max(dists)) #bw="ucv" --> tres lent !

	upD = min(k, n-1, ncol(data))
    for (d in 2:upD) {
		svols = sampleVols(N, d, densit)
		cvols = c()
        for (i in 1:N) {
            #pick d+1 points at random in a neighborhood
            p = sample(1:n, 1)
			inds = sample(1:(k+1), d+1)
            locData = data[ NI[[p]] [inds], ]
            cvols = c(cvols, simpVol(locData))
        }

		m = min(svols,cvols)
		M = max(svols,cvols)
		d1 = density(svols, kernel = "epanechnikov", from=m, to=M)$y
		d2 = density(cvols, kernel = "epanechnikov", from=m, to=M)$y
		delta = sqrt( (1/length(d1)) * sum( (d2-d1)^2 ) )
		if (delta > thtest) return (d-1)
    }
	
    return (upD)
}
