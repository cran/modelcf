#main function RML (Riemnannian Manifold Learning)
RML = function(data, d, kmin=0, kmax=0, tsoft=0.1, thlvl=0.3, hdth=2, advhr=FALSE) {
    n = nrow(data)

    #step1 : neighborhoods
    # get knn(s) :
    if (kmax==0) kmax = max( min(n-1, 2*getKnn(n)), adaptAux(data,d,0.95)$kmax )
    if (kmin==0) kmin = 1
    else kmax = max(kmax, kmin)
	distsEuc = as.matrix(dist(data)) ; diag(distsEuc) = Inf

	#compute rough graph distances (if advanced heuristic asked) :
	rgdists = -1.0
	
	if (advhr) {
		rgdists = matrix(0.0,nrow=n,ncol=n)
		NI_dij = simpleNeighbs(data, getLimitConnex(data) )
		rgLocDists = matrix(0.0,nrow=n,ncol=n)
		#keep only graph neighborhoods :
		for (i in 1:n) {
			rgLocDists[ i,NI_dij[[i]] ] = distsEuc[ i,NI_dij[[i]] ]
		}
		rgLocDists = pmax(rgLocDists, t(rgLocDists))
		for (i in 1:n) {
			dj_tmp = .C("dijkstra",dists=as.double(rgLocDists),n=as.integer(n),start=as.integer(i),
							geodsPredsLevels=double(3*n), PACKAGE=pkgnm())$geodsPredsLevels
			rgdists[i,] = dj_tmp[1:n]
		}
	}

    #get neighborhoods :
    NI = neighbs_RML(data, rgdists, kmin, kmax, tsoft)
    NI = testConnexity(data,NI,kmin)
	
    ######################################################
    #step2 : Dijkstra from the "average" function
    locDists = matrix(0.0,nrow=n,ncol=n)
    #keep only graph neighborhoods :
    for (i in 1:n) {
        locDists[ i,NI[[i]] ] = distsEuc[ i,NI[[i]] ]
    }
    locDists = pmax(locDists, t(locDists))

    #compute (approximate) geodesic distances from each point :    
    bestDijkVect = 0 ; bestInd = 1 ; bestRad = Inf
    seqVect = 1:n
    for (i in 1:n) {
        dj_tmp = .C("dijkstra",dists=as.double(locDists),n=as.integer(n),start=as.integer(i),
                    geodsPredsLevels=double(3*n), PACKAGE=pkgnm())$geodsPredsLevels
        s = max(dj_tmp[1:n])
        if (s < bestRad) {
            bestRad = s
            bestDijkVect = dj_tmp
            bestInd = i
        }
    }
    med = bestInd #origin index of the smallest radius

    #retrieve predecessors and levels :
    preds = bestDijkVect[(n+1):(2*n)]
    levels = bestDijkVect[(2*n+1):(3*n)]
    tre = hdth
    if (tre == 0) tre = max(1, floor(thlvl * max(levels)) )

    #preliminary step : compute base vectors of the tangent space at data[med,], in columns
    medNeighbs = t( t(data[ NI[[med]], ]) - data[med,])
    s = svd(medNeighbs,nu=0,nv=d)
    basis = t(as.matrix(s$v))


    ######################################################
    #step3 : find the global coordinates for all points
    embedding = matrix(nrow=n,ncol=d) ; embedding[med,] = rep(0.0,d)
    coordsOK = rep(FALSE,n) ; coordsOK[med] = TRUE

	homt_vect = c()
	orig_embs_inds = c()
	
    #easy case : neighbor of the starting point
    for (i in 1:n) {
        if (levels[i] > 0 && levels[i] <= tre) {
            X = as.double( basis %*% (data[i,] - data[med,]) ) 
			orig_embs_inds = c(orig_embs_inds,i)
			homt = sqrt(sum(X^2)) / sqrt(sum((data[i,]-data[med,])^2))
			homt_vect = c(homt_vect, homt)
            embedding[i,] = (1.0 / homt) * X
            coordsOK[i] = TRUE
        }
    }
	
	#locally learn unscaled embeddings ([generalized] linear regression)
	x = embedding[orig_embs_inds,]
	x = cbind( rep(1,length(orig_embs_inds)), x, x^2 )
	y = homt_vect
	invtxx = mppsinv( t(x) %*% x )
	linCoefs = invtxx %*% t(x) %*% y

    #then, while all coordinates haven't been found..
    currentLevel = tre + 1
    avknn = round( 0.5* (kmin+kmax) )
    while (sum(coordsOK) < n) {
        for (i in 1:n) {
            #second case : we get the embedding by numerical optimization (differs from article)
            if (levels[i] == currentLevel) {
                dists = distsEuc[preds[i],coordsOK==TRUE] #distances from predecessor to already treated points
                ixDone = seqVect[coordsOK==TRUE] #corresponding indices
                nbContraintes = min(avknn, length(ixDone)-1) #reasonable number of constraints (?!)
                srt = sort(dists,index.return=TRUE)
                indsN = srt$ix[1:nbContraintes] #"normal" indices corresponding to the nbContraintes nearest
                indices = ixDone[indsN] #back to true indices

                #call numerical "solver" :
                embedding[i,] = solveConstraints(data, embedding, i, preds[i], indices, d)
                coordsOK[i] = TRUE
            }
        }
        currentLevel = currentLevel + 1
    }

    return (list("data"=data,"embed"=embedding,"lincofs"=linCoefs,"lvls"=levels,"fdata"=data[med,],"fbasis"=basis,"tre"=tre,"NI"=NI))
}


########################################################
## RML reconstruction :

#warning : n > 4 (should always be the case)
RML_rec = function(RLout,newEmb) {
    data = RLout$data
    n = nrow(data)

    # STEP 0: COMPUTE PAIRWISE DISTANCES and FIND NEIGHBORS
    embeds = as.matrix(RLout$embed)
	d = ncol(embeds)
    distsN = sqrt( colSums( (t(embeds) - newEmb)^2 ) )
    srtN = sort(distsN,index.return=TRUE)

    if (srtN$x[1] < EPS() ) return ( data[srtN$ix[1],] ) #shortcut
    closest = srtN$ix[1]
    rlK = max( 2, d, length(RLout$NI[[closest]]) )

    #case 1 : can be approximated by local (origin) basis
    if (sum(RLout$lvls[ srtN$ix[1:rlK] ] <= RLout$tre) >= ceiling(rlK/1.5)) {
        #neighbors of level <=tre are the most numerous :
        tmp = RLout$fdata + ( sum( c(1, newEmb, newEmb^2) * RLout$lincofs) * newEmb ) %*% RLout$fbasis
        return (as.double(tmp))
    }

    #retrieve kNN neighborhood "in RML order" from z_p (if enough..)
    dists = as.matrix(dist(embeds))
    diag(dists) = Inf
    srt = sort(dists[closest,],index.return=TRUE)
    rlNeighbs = RLout$NI[[closest]]
    if (length(rlNeighbs) < rlK) {
		isIn = rep(FALSE,n)
		isIn[rlNeighbs] = TRUE
		for (i in 1:(n-1)) {
			if ( ! isIn[srt$ix[i]] ) rlNeighbs = c(rlNeighbs, srt$ix[i])
			if (length(rlNeighbs) == rlK) break
		}
	}
    
    #case level > tre : we start with a local PCA (knn for ACP computation = ?)
    #the nearest neighbor is used instead of the mean :
    ctrData = t( t(data[srt$ix[1:(rlK+1)],]) - data[closest,])
    s = svd(ctrData,nu=0,nv=d) ; V = as.matrix(s$v)
    locData = t( t(data[rlNeighbs,]) - data[closest,])
    redCoords = locData %*% V

    #call numerical "solver" :
    data_ = rbind( embeds[closest,], as.matrix(embeds[rlNeighbs,]), newEmb )
    embedding_ = rbind(rep(0.0,d), redCoords, rep(0.0,d))
    rlK = length(rlNeighbs)
    tmp = solveConstraints(data_,embedding_,rlK+2,1,2:(rlK+1),d)
    res = data[closest,] + V %*% tmp #preserve distances
    return (res)
}

