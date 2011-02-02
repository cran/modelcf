#function LPcaML : local PCA + alpha-TSLN
LPcaML = function(data, d, adn="none", k=0, alpha=0.5, trcv=0.7) {
	require(class)
    n=nrow(data) ; D = ncol(data)
    minSizeV = min(n, ceiling((1.0/alpha)*(d+1)) )
    if (k==0) k = getKnn(n)
    k_orig = k
    k = max(minSizeV-1,k)

    #ensure connexity :
    NI = getNI(data, adn, d, k)
    NI = testConnexity(data, NI, k)
    NI = getMutual(NI, TRUE)
    # add "self" indices :
    for (i in 1:n) NI[[i]] = c(i,NI[[i]])
    ks = sapply(NI, length)

    #initialisation : seek for the point where manifold is the most dense
    #we could try to find the "most linear" point, but this is very costly..
    bestInd = 0
    minDist = Inf
    for (i in 1:n) {
        tmp = sqrt( sum (data[i,] - data[ NI[[i]] [ks[i]], ])^2 )
        if (tmp < minDist) {
            minDist = tmp
            bestInd = i
        }
    }

    #step 1 : finding an alpha-TSLN from bestInd.
    #loop over all neighborhoods and include them if they overlap ;
    sets = rep(0,n) #indices of sets of indices already treated, wrt. alpha-TSLN order
    sets[1] = bestInd
    current = rep(FALSE,n) #the indices of current sets
    current[ NI[[bestInd]] ] = TRUE
    curIndSet=1
    
    while (sum(current) < n) {    
        bestCand = 0 ; bestValCand = Inf
        for (i in 1:n) {
            tmp = sum(current[ NI[[i]] ])
            deltaAbs = abs(tmp - alpha*ks[i])
            if (tmp < ks[i] && deltaAbs < bestValCand) {
                bestCand = i
                bestValCand = deltaAbs
            }
        }
        current[ NI[[bestCand]] ] = TRUE
        curIndSet = curIndSet + 1
        sets[curIndSet] = bestCand
    }
    sets = sets[sets > 0]
    nbSets = length(sets)

    #step 2 : compute local PCA's to get all local coefficients +
    #find affine transformations to get global coords.
    locCoords = matrix(0.0,nrow=n,ncol=d)
    globCoords = matrix(0.0,nrow=n,ncol=d)
    indSets = as.list(rep(0,nbSets)) #sets of indices of "alpha bubbles" (in order)
    indSets[[1]] = NI[[bestInd]]
    firstInds = NI[[bestInd]]
    means = as.list(rep(0,nbSets)) #means of "alpha-bubbles"
    bases = as.list(rep(0,nbSets)) #corresponding local basis
    invB = as.list(rep(0,nbSets)) #storage of inverse B matrices
    invB[[1]] = diag(ks[bestInd]) #first one is identity (and useless)

    means[[1]] = colMeans(data[firstInds,])
    ctrData = t( t(data[firstInds,]) - means[[1]] )
    s = svd(ctrData,nu=0,nv=d)
    bases[[1]] = t(as.matrix(s$v))
    locCoords[firstInds,] = ctrData %*% s$v
    globCoords[firstInds,] = locCoords[firstInds,] #anchor points
    
    #loops to get all coordinates patches
    for (i in 2:nbSets) {

        if (i > nbSets) break; #bad : only one set we shouldn't loop !

        curInds = NI[[ sets[i] ]]
        indSets[[i]] = NI[[ sets[i] ]]

        #assignment of local coordinates :
        means[[i]] = colMeans(data[curInds,])
        ctrData = t( t(data[curInds,]) - means[[i]] )
        s = svd(ctrData,nu=0,nv=d)
        locCoords[curInds,] = ctrData %*% s$v
        bases[[i]] = t(as.matrix(s$v))

        #find the union of C_i inter C_j, j<i :
        inters = rep(FALSE,n)
        for (j in 1:(i-1)) inters[ intersect(curInds,NI[[ sets[j] ]]) ] = TRUE

        #find the matrix B transforming local into global coordinates
        # ==> standard linear regression
        sint = sum(inters)
        if (sint>1) unYp = cbind(rep(1,sint),locCoords[inters,])
        else unYp = t(as.matrix(c(1,locCoords[inters,])))

        #solve (1 Y')B=Y (see article)
        B = matrix()
        if (sint >= d+1) B = mppsinv(unYp) %*% globCoords[inters,] #always this case.. (in principle)

		else if (sint == 1) B = rbind(globCoords[inters,],matrix(0,nrow=d,ncol=d)) #shouldn't happen, but...

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
        else invB[[i]] = rbind(B[1,],matrix(0.0,nrow=d,ncol=d))
    }

	#final stage : train a classifier z_i --> labels
	labels = rep(0,n)
	for (i in nbSets:1) labels[ NI[[ sets[i] ]] ] = i
	labels = reordering(labels)

    #optimize model parameters (number of neighbors) :
    kclass = optimParams_classif(globCoords,labels,"kNN",floor(trcv*n),trcv)
    #build model :
    classifObj = learnClassif(globCoords, labels, "kNN", kclass)

    return (list("data"=data,"embed"=globCoords,"invB"=invB,"indSets"=indSets,"means"=means,"bases"=bases, "classi"=classifObj))
}


#################################
## adapted reconstruction :

#output a curve from some d-dimensional representation
LPcaML_rec = function(LPout,newEmb) {
	d = ncol(LPout$embed)
    #predict classification (=="alpha-bubble") :
	cl = predictClassif(LPout$classi,newEmb)

	if (cl>1) locCoords = (newEmb-LPout$invB[[cl]][1,]) %*% LPout$invB[[cl]][2:(d+1),]
    else locCoords = newEmb
    reconst = LPout$means[[cl]] + linear_rec(LPout$bases[[cl]], locCoords)
    return (reconst)
}

