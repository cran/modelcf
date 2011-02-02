## clustering of inputs (or-) outputs data

#wrapper for hierarchical clustering, Ward linkage :
phclust = function(dissims, K) {
    hct = hclust(as.dist(dissims),method="ward")
    return ( cutree(hct, K) )
}

#main function :
gtclusts = function(method, data, K, d=min(5,ncol(data)), adn="none", knn=0, symm=TRUE, weight=FALSE, sigmo=FALSE) {

	checkDependP(mclust=method)
    if (knn==0) knn = getKnn(nrow(data))

    dissims = matrix()
    if (method=="HDC") dissims = hitorct(data, adn, d, knn, FALSE, symm, weight)
    else if (method=="CTH" || method=="CTKM") dissims = ctdists(data, adn, d, knn, symm, weight, sigmo)
    else if (method=="CH") dissims = as.matrix(dist(data))

    clusts = c()
    if (method=="HDC" || method=="CTKM") clusts = km_dists(dissims, K)
    else if (method=="CTH" || method=="CH") clusts = phclust(as.dist(dissims), K)

    else if (method=="PCA") {
        #always the simplified version (gives good results); the other one is too costly..
        clusts = km_PCA(data, K, d, TRUE)
    }
    else if (method=="KM") clusts = (kmeans(data,K,iter.max=100,nstart=10))$cluster
	else if (method=="spec") clusts = as.integer(specc(data,K,kernel="rbfdot",kpar="automatic"))

    return (clusts)
}

#procedure to determine K by stability and prediction, and return the corresponding partition :
findK_gtclusts = function(x, y, method, d=min(10,ncol(x)), adn="none", knn=0,symm=TRUE,weight=FALSE,sigmo=FALSE,
                          minszcl=30,maxcl=Inf,mclass="kNN",taus=0.8,Ns=10,tauc=0.8,Nc=10,trcv=0.7,nstagn=10) {

	checkDependP(mclust=method,mclass=mclass)
    n = nrow(y) ; seqVect = 1:n
    lte = n - min( n-2, max(2,round(trcv*n)) )
    ltr = n - lte
    if (knn==0) knn = getKnn(n)

    maxcl = min(n-1,maxcl)
    knnTR = getKnn(ltr)

    oldPref = rep(1, n)
    lastK = 1
    K = 2
    while (K <= maxcl) {
        P_ref = reordering( fusion_smcl( y, gtclusts(method, y, K, d, adn, knn, symm, weight, sigmo), minszcl ) )

        if (max(P_ref) > max(oldPref)) {
            lastK = K
            tau = 0.0

            #check stability to sub-sampling :
            for (i in 1:Ns) {
                linter = 0
                while (linter == 0) {
                    U1 = sample(seqVect,ltr)
                    U2 = sample(seqVect,ltr)
                    I = intersect(U1, U2)
                    linter = length(I)
                }
                knnI = getKnn(linter)
                indsU1 = rep(FALSE, n)
                indsU1[U1] = TRUE
                indsU1[I] = FALSE
                indsU2 = rep(FALSE, n)
                indsU2[U2] = TRUE
                indsU2[I] = FALSE
                reorderY1 = rbind(y[I,],y[indsU1,])
                reorderY2 = rbind(y[I,],y[indsU2,])
				locMaxszcl = ceiling(minszcl*length(I)/n)
                P1 = reordering( fusion_smcl( reorderY1, gtclusts(method, reorderY1, K, d, adn, knnI, symm, weight, sigmo), locMaxszcl ) )
                P2 = reordering( fusion_smcl( reorderY2, gtclusts(method, reorderY2, K, d, adn, knnI, symm, weight, sigmo), locMaxszcl ) )
                set1 = P1[1:linter]
                set2 = P2[1:linter]

                #we compute the rate of "same cluster" elements with Meila (2002) variation of information index
                tau = tau + varInfo(set1,set2) / Ns
            }

            if (tau < taus) return (oldPref)
            tau = 0.0

            #cross-validation (prediction) ; measures also the stability to subsampling
            for (i in 1:Nc) {

                #cluster a training set :
                testInds = sample(seqVect,lte)
                trainInds = seqVect[-testInds]
                P1 = fusion_smcl( y[trainInds,], gtclusts(method, y[trainInds,], K, d, adn, knnTR, symm, weight, sigmo), ceiling(minszcl*trcv) )

                #optimize model parameters :
                opt = optimParams_classif(x[trainInds,],P1,mclass,knnTR,trcv)
                #build model :
                classifObj = learnClassif(x[trainInds,], P1, mclass, opt)
                #predict outputs classification :
                clElems = predictClassif(classifObj,as.matrix(x[testInds,]))

                #now compare it to the classification in P2 (==P_ref) :
                set1 = reordering(clElems)
                set2 = P_ref[testInds]                
                maxInd1 = max(set1)
                maxInd2 = max(set2)

                #renumbering step to match set1 clusters to set2 the best possible way, using hungarian algorithm
                #"utils" are simply the cardinals of clusters intersections
                #WARNING : set1 can have less clusters than set2, but not more
                assign_cl = c()
                if (maxInd1 <= maxInd2) {    
                    assign_cl = .C("getOptAssign",P1=as.integer(set1),P2=as.integer(set2),maxInd=as.integer(maxInd2),
                                   n=as.integer(lte),maxim=as.integer(1),assign_cl=integer(maxInd2),PACKAGE=pkgnm())$assign_cl
                    Ptmp = set1
                    for (i in 1:maxInd1) set1[Ptmp == i] = assign_cl[i]+1
                }
                else {
                    assign_cl = .C("getOptAssign",P1=as.integer(set2),P2=as.integer(set1),maxInd=as.integer(maxInd1),
                                   n=as.integer(lte),maxim=as.integer(1),assign_cl=integer(maxInd1),PACKAGE=pkgnm())$assign_cl
                    Ptmp = set2
                    for (i in 1:maxInd2) set2[Ptmp == i] = assign_cl[i]+1
                }

                #we compute the rate of well-classified elements
                tau = tau + sum(set1 == set2) / (Nc*lte)
            }

            if (tau < tauc) return (oldPref)
            oldPref = P_ref
        }
        K = K+1
        if (K > lastK+nstagn || K >= ltr) return (oldPref)
    }

    return (oldPref)
}

################################################################
# main method : call the previous one on outputs, and then on each inputs cluster.
gtclusts_inout = function(x, y, method, d=min(10,ncol(x)), redy=FALSE, adn="none", knn=0, symm=TRUE, weight=FALSE,sigmo=FALSE,
                          minszcl=30,maxcl=Inf, mclass="kNN", taus=0.8, Ns=10, tauc=0.8, Nc=10, trcv=0.7, verb=TRUE, nstagn=10) {

	checkDependP(mclust=method,mclass=mclass)
	if (d==0) d = dest_levi(y, knn)
    if (redy) {    
        #basic outputs transformation with PCA to fasten overall process :
        s = svd( y, nu=0, nv= min(2*d, ncol(y)) )
        y = y %*% s$v
    }

    if (verb) cat("    --- outputs ---\n")
    gcOut = findK_gtclusts(x, y, method, d, adn, knn, symm, weight, sigmo, minszcl, maxcl, mclass, taus, Ns, tauc, Nc, trcv, nstagn)
    nbC_out = max(gcOut)
    if (nbC_out==1 || nbC_out == maxcl) return(gcOut)

    finalClusts = gcOut ; nbClusts = nbC_out
    if (verb) cat("    --- inputs ---\n")
    for (i in 1:nbC_out) {
        inaout = as.matrix(x[gcOut==i,])
        gcIn = findK_gtclusts(inaout,inaout,method,d,"none",0,symm,weight,sigmo,minszcl,maxcl-nbClusts+1,mclass,taus,Ns,tauc,Nc,trcv,nstagn)
        finalClusts[gcOut==i] = max(finalClusts) + gcIn
        nbClusts = nbClusts + max(gcIn) - 1
    }
    return ( reordering(finalClusts) )
}

