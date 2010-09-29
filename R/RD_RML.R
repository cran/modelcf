#main function RML (Riemnannian Manifold Learning)
RML = function(data, d, adn=FALSE, knnmin=0, knnmax=0, tsoft=0.1, thlvl=0.3, hdth=0) {
    n = nrow(data)

    #step1 : neighborhoods
    # get knn(s) :
    if (knnmin==0) knnmin = max(4, d+1)
    if (knnmax==0) {
        if (adn) knnmax = adaptiveNeighbs1(data,d)
        else knnmax = getKnn(n)
    }
    knnmax = max(knnmin,knnmax)

    #get neighborhoods :
    NI = neighbs_RML(data, knnmin, knnmax, tsoft)
    NI = testConnexity(data,NI,knnmin)


    ######################################################
    #step2 : Dijkstra from the "average" function
    distsEuc = as.matrix(dist(data)) ; diag(distsEuc) = Inf
    locDists = matrix(0.0,nrow=n,ncol=n)
    #keep only graph neighborhoods :
    for (i in 1:n) {
        locDists[ i,NI[[i]] ] = distsEuc[ i,NI[[i]] ]
    }
    locDists = pmax(locDists, t(locDists))

    #compute (approximate) geodesic distances from each point :
    bestDijkVect = 0 ; bestInd = 1 ; bestSumDists = Inf
    seqVect = 1:n
    for (i in 1:n) {
        dj_tmp = .C("dijkstra",dists=as.double(locDists),n=as.integer(n),start=as.integer(i),
                    geodsPredsLevels=double(3*n), PACKAGE=pkgnm())$geodsPredsLevels
        s = sum(dj_tmp[1:n])
        if (s < bestSumDists) {
            bestSumDists = s
            bestDijkVect = dj_tmp
            bestInd = i
        }
    }
    med = bestInd #origin index of the smallest sum of distances

    #retrieve predecessors and levels :
    preds = bestDijkVect[(n+1):(2*n)]
    levels = bestDijkVect[(2*n+1):(3*n)]
    tre = hdth
    if (tre == 0) tre = max(1, ceiling(thlvl * max(levels)) )

    #preliminary step : compute base vectors of the tangent space at data[med,], in columns
    medNeighbs = t( t(data[ NI[[med]], ]) - data[med,])
    s = svd(medNeighbs,nu=0,nv=d)
    basis = as.matrix(s$v)


    ######################################################
    #step3 : find the global coordinates for all points
    embedding = matrix(nrow=n,ncol=d) ; embedding[med,] = rep(0.0,d)
    coordsOK = rep(FALSE,n) ; coordsOK[med] = TRUE

    #easy case : neighbor of the starting point
    for (i in 1:n) {
        if (levels[i] > 0 && levels[i] <= tre) {
            X = qr.solve(basis,data[i,]-data[med,])
            embedding[i,] = ( sqrt(sum((data[i,]-data[med,])^2)) / sqrt(sum(X^2)) ) * X
            coordsOK[i] = TRUE
        }
    }

    #then, while all coordinates haven't been found..
    currentLevel = tre + 1
    avknn = round( 0.5* (knnmin+knnmax) )
    while (sum(coordsOK) < n) {
        for (i in 1:n) {
            #second case : we get the embedding by numerical optimization (differs from article)
            if (levels[i] == currentLevel) {
                dists = distsEuc[preds[i],coordsOK==TRUE] #distances from predecessor to already treated points
                ixDone = seqVect[coordsOK==TRUE] #corresponding indices
                nbContraintes = min(avknn, length(ixDone)-1) #reasonable number of constraints (?!)
                srt = sort(dists,index.return=TRUE)
                indsN = srt$ix[2:(nbContraintes+1)] #"normal" indices corresponding to the nbContraintes nearest
                indices = ixDone[indsN] #back to true indices
                indices = indices[indices != med] #..WITHOUT med if present (otherwise A_rmlaux will contain NaN)

                #call numerical "solver" :
                embedding[i,] = solveConstraints(data, embedding, i, preds[i], indices, d)
                coordsOK[i] = TRUE
            }
        }
        currentLevel = currentLevel + 1
    }

    return (list("data"=data,"embed"=embedding,"lvls"=levels,"fdata"=data[med,],"fbasis"=basis,"tre"=tre,"tsoft"=tsoft,"NI"=NI))
}


########################################################
## RML reconstruction :

#warning : n > 4 (should be always the case)
RML_rec = function(RLout,newEmb) {

    data = RLout$data
    d = ncol(RLout$embed)
    n = nrow(data)

    # STEP 0: COMPUTE PAIRWISE DISTANCES and FIND NEIGHBORS
    embeds = RLout$embed
    dists = sqrt( colSums( (t(embeds) - newEmb)^2 ) )
    srt = sort(dists,index.return=TRUE)

    if (srt$x[1] < EPS() ) return ( data[srt$ix[1],] ) #shortcut
    #knn is taken as the number of neighbors of the nearest embedded points
    knn = length(RLout$NI[[ srt$ix[1] ]])
    neighbs = srt$ix[1:knn]

    #compute kNN neighborhood "in RML order" [could be useful if extended with less constraints..]
    rlNeighbs = c(neighbs[1])
    addIndBool = rep(FALSE,knn) ; addIndBool[1] = TRUE
    for (i in 2:knn) {
        #computation of dot products
        embrl = as.matrix(embeds[rlNeighbs,])
        cs = colSums( (t(embrl) - newEmb) * (t(embrl) - embeds[neighbs[i],]) )
        if (min(cs) >= RLout$tsoft) {
            rlNeighbs = c(rlNeighbs,neighbs[i])
            addIndBool[i] = TRUE
        }
    }
    for (i in 2:knn) {
        if (!addIndBool[i]) rlNeighbs = c(rlNeighbs,neighbs[i])
    }

    seqVect = 1:n
    #case 1 : can be approximated by local (origin) basis
    if (sum(RLout$lvls[ seqVect[rlNeighbs] ] <= RLout$tre) >= ceiling(knn/2.0)) {
        #neighbors of level <=tre are the most numerous :
        tmp = RLout$fdata + newEmb %*% t(RLout$fbasis)
        return (as.double(tmp))
    }

    #case level > tre : we start with a local PCA
    #trick : the nearest neighbor is used instead of the mean :
    ctrData = t( t(data[rlNeighbs,]) - data[rlNeighbs[1],])
    s = svd(ctrData,nu=0,nv=d) ; V = as.matrix(s$v)
    redCoords = ctrData %*% V
    #call numerical "solver" :
    data_ = rbind( as.matrix(embeds[rlNeighbs,]), newEmb )
    embedding_ = rbind(redCoords,rep(0.0,d))

    tmp = solveConstraints(data_,embedding_,knn+1,1,2:(knn-1),d)
    res = data[rlNeighbs[1],] + V %*% tmp #preserve distances
    return (res)
}

