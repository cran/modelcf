#~fusion RML + LPcaML
GCEM = function(data, d, adn=FALSE, knn=0, alpha=0.5, tsoft=0.1, thlvl=0.3, hdth=0) {
    n=nrow(data) ; D = ncol(data)
    minSizeV = min(n, ceiling((1.0/alpha)*(d+1)) )
    if (knn==0) knn = getKnn(n)
	kNN_orig = knn
    knn = max(minSizeV-1,knn)

    #ensure connexity :
    NI = getNI(data, adn, d, knn)
    NI = testConnexity(data, NI, knn)
    # add "self" indices :
    for (i in 1:n) NI[[i]] = c(i,NI[[i]])
    knns = sapply(NI, length)

    #initialisation : seek for the point where manifold is the most dense
    #we could try to find the "most linear" point, but this is very costly..
    bestInd = 0
    minDist = Inf
    for (i in 1:n) {
        tmp = sqrt( sum (data[i,] - data[ NI[[i]] [knns[i]], ])^2 )
        if (tmp < minDist) {
            minDist = tmp
            bestInd = i
        }
    }

    #step 1 : finding an alpha-TSLN from bestInd.
    #loop over all neighborhoods and include them if they overlap ;
    sets = rep(0,n) #indices of sets of indices already treated, wrt. alpha-TSLN order
    sets[1] = bestInd
    current = rep(1,n) #the NON-indices of current sets (1 if NOT current)
    current[ NI[[bestInd]] ] = 0
    curIndSet=1

    while (sum(current) > 0) {
        bestCand = 0 ; bestSumCand = 0
        savCand = 0 ; bestSav = Inf
        for (i in 1:n) {
            tmp = sum(current[ NI[[i]] ])
            if (tmp > bestSumCand) {
                if (tmp <= (1.0-alpha)*knns[i]) {
                    bestCand = i
                    bestSumCand = tmp
                }
                else {
                    if (tmp < bestSav) {
                        savCand = i
                        bestSav = tmp
                    }
                }
            }
        }
        if (bestCand==0) bestCand = savCand
        current[ NI[[bestCand]] ] = 0
        curIndSet = curIndSet + 1
        sets[curIndSet] = bestCand
    }
    sets = sets[sets > 0]
    nbSets = length(sets)


    #step 2 : compute local RML's to get all local coefficients +
    #find affine transformations to get global coords.
    locCoords = matrix(0.0,nrow=n,ncol=d)
    globCoords = matrix(0.0,nrow=n,ncol=d)
    indSets = as.list(rep(0,nbSets)) #sets of indices of "alpha bubbles" (in order)
    indSets[[1]] = NI[[bestInd]]
    firstInds = NI[[bestInd]]
    invB = as.list(rep(0,nbSets)) #storage of inverse B matrices
    invB[[1]] = diag(knns[bestInd]) #first one is identity (and useless)

    locKnn = min( knns[bestInd]-1, max( 1, min(30,knn%/%3) ) )
    rlObjects = as.list(rep(0,nbSets))
    rl = RML(data[firstInds,],d,FALSE,tsoft,thlvl,hdth,locKnn,locKnn)
    rlObjects[[1]] = rl
    locCoords[firstInds,] = rl$embed
    globCoords[firstInds,] = locCoords[firstInds,] #anchor points

    #loops to get all coordinates patches
    for (i in 2:nbSets) {

        if (i > nbSets) break; #bad : only one set we shouldn't loop !

        curInds = NI[[ sets[i] ]]
        indSets[[i]] = NI[[ sets[i] ]]

        #assignment of local coordinates :
		locKnn = min( knns[sets[i]]-1, max( 1, min(30,knn%/%3) ) )
        rl = RML(data[curInds,],d,FALSE,tsoft,thlvl,hdth,locKnn,locKnn)
        rlObjects[[i]] = rl
        locCoords[curInds,] = rl$embed

        #find the union of C_i \inter C_j, j<i :
        inters = rep(FALSE,n)
        for (j in 1:(i-1)) inters[ intersect(curInds,NI[[ sets[j] ]]) ] = TRUE

        #find the matrix B transforming local into global coordinates
        # ==> standard linear regression
        B = matrix()
        sint = sum(inters)
        if (sint>1) unYp = cbind(rep(1,sint),locCoords[inters,])
        else unYp = t(as.matrix(c(1,locCoords[inters,])))

        #solve (1 Y')B=Y (see article)
        if (sint >= d+1) B = mppsinv(unYp) %*% globCoords[inters,]

        #a little more complicated :
        else {
            unYp_square = unYp[,1:nrow(unYp)]
            B = mppsinv(unYp_square) %*% globCoords[inters,]
            B = rbind(B,matrix(0,nrow=d+1-nrow(unYp),ncol=d))
        }

        #apply B to get global coordinates :
        allF = rep(FALSE,n)
        allF[curInds] = TRUE ; allF[inters] = FALSE #symmetric difference
        seqVect = 1:n
        estimGlob = seqVect[allF==TRUE] #indices to be estimated
        if (length(estimGlob) > 1)
            globCoords[estimGlob,] = cbind(rep(1,length(estimGlob)),locCoords[estimGlob,]) %*% B
        else globCoords[estimGlob,] = c(1,locCoords[estimGlob,]) %*% B

        #find B pseudo-inverse (for reconstruction) :
        s = svd(B[2:(d+1),]) ; tmpD=s$d
        cutDiag = tmpD[tmpD > EPS()]
        if (length(cutDiag) > 0) {
            invD = matrix(0.0,nrow=d,ncol=d)
            diag(invD)[1:length(cutDiag)] = 1.0 / cutDiag
            invB[[i]] = rbind(B[1,], s$v %*% invD %*% t(s$u))
        }
        else invB[[i]] = rbind(B[1,],matrix(0,nrow=d,ncol=d))
    }

    return (list("data"=data,"embed"=globCoords,"invB"=invB,"indSets"=indSets,"rlObjs"=rlObjects,"knn"=kNN_orig))
}


########################################################
## GCEM reconstruction :

#output a curve from some d-dimensional representation
GCEM_rec = function(GCout,newEmb) {
    data = GCout$data
    n = nrow(data) ; D = ncol(data)
    knn = GCout$knn
    if (knn==0) knn = getKnn(n)
    d = ncol(GCout$embed)

    # STEP 0: COMPUTE PAIRWISE DISTANCES and FIND NEIGHBORS
    embeds = GCout$embed
    dists = sqrt( colSums( (t(embeds) - newEmb)^2 ) )
    srt = sort(dists,index.return=TRUE)

    if (srt$x[1] < EPS() ) return ( data[srt$ix[1],] ) #shortcut
    neighbs = srt$ix[1:knn]

    #histogram of sets indices among neighbors :
    indSets = GCout$indSets ; nbSets = length(indSets)
    countSets = rep(0,nbSets)
    for (i in 1:knn) {
        for (j in 1:nbSets) {
            if (length(intersect(neighbs[i],indSets[[j]])) > 0)
                countSets[j] = countSets[j] + 1
        }
    }
    reconst = rep(0.0,D)
    wm = which.max(countSets) ; seqVect = 1:nbSets
    exaequos = seqVect[ countSets == countSets[wm] ] #indices of "winning" neighbors

    for (i in exaequos) {
        if (i>1) locCoords = (newEmb-GCout$invB[[i]][1,]) %*% GCout$invB[[i]][2:(d+1),]
        else locCoords = newEmb
        tmp = RML_rec(GCout$rlObjs[[i]],as.double(locCoords))
        reconst = reconst + tmp
    }

    reconst = reconst / length(exaequos)
    return (reconst)
}

