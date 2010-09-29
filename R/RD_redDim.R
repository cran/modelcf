##################################################################
## calls of functions for dimension reduction and reconstruction

#dimensionality reduction
nlin_redDim = function(method, data, d, adn, knn, alpha, knnmin, knnmax, tsoft, thlvl, hdth) {

    if (knn==0) knn = getKnn(nrow(data))
    if (method=="LPcaML") return ( LPcaML(data, d, adn, knn, alpha) )
    if (method=="RML") return ( RML(data,d,adn,knnmin,knnmax,tsoft,thlvl,hdth) )
    if (method=="GCEM") return ( GCEM(data, d, adn, knn, alpha, tsoft, thlvl, hdth) )
}

#reconstruction
nlin_adaptRec = function(method, embobj, newEmb) {
    fRec = match.fun(paste(method,"_rec",sep=""))
    return ( fRec(embobj,newEmb) )
}

