########################################################
## two methods for dimension estimation :

#return local estimated dimension
locdim1 = function(x, knn) {
    refVect = 1:length(x)
    classI = rank(x,ties.method="random")
    j1 = refVect[classI == knn%/%2]
    j2 = refVect[classI == knn]
    return ( log(2) / (log(x[j2] / x[j1])) )
}

#dimension estimation according to Farahmand et al. (2006)
dimest1 = function(data, knn) {
    dists = as.matrix(dist(data))
    diag(dists) = Inf
    knn = max(2,knn)
    dim = 0.0

    #compute all local dimensions and average them
    dim = mean( apply(dists, 1, locdim1, knn) )
    #return the nearest integer
    return ( round(dim) )
}

#simplex approximation of the dimension (Lin and Zha, RML 2006)
dimest2 = function(data, knnmin, knnmax, tsoft=0.1) {
    n = nrow(data)
    prodTab = n*knnmax

    # initial step : compute the local neighborhoods [semi-adaptative, like in RML]
    tabK = .C("estimK_RML",data=as.double(t(data)),n=as.integer(n),m=as.integer(ncol(data)),
            knnmin=as.integer(knnmin),knnmax=as.integer(knnmax),tsoft=as.double(tsoft),
            tabK=integer(prodTab), PACKAGE=pkgnm())$tabK

    #call C routine to end the job :
    dim = .C("simpDim",neighbs=as.integer(tabK),n=as.integer(n),kNN=as.integer(knnmax),
            dim=integer(1), PACKAGE=pkgnm())$dim
    return (dim)
}

