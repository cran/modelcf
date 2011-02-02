##################################################################
## calls of functions for dimension reduction and reconstruction

#dimensionality reduction
nlin_redDim = function(method, data, d, adn, k, alpha, trcv, kmin, kmax, tsoft, thlvl, hdth, advhr) {

    if (k==0) k = getKnn(nrow(data))

    if (method=="LPcaML") return ( LPcaML(data, d, adn, k, alpha, trcv) )
    if (method=="RML") return ( RML(data, d, kmin, kmax, tsoft, thlvl, hdth, advhr) )
}

#reconstruction
nlin_adaptRec = function(method, embobj, newEmb) {
    fRec = match.fun(paste(method,"_rec",sep=""))
    return ( fRec(embobj,newEmb) )
}

