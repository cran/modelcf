##############################################
## ("linear") orthonormal basis

# select the sub-basis of most variable coefficients
basisMostVar = function(basis,data,nbCoefs) {
    F = data %*% t(basis)
    rnCoefs = ncol(F)

    if (nbCoefs >= rnCoefs) return (basis)
    V = rep(0,rnCoefs)
    for (i in 1:rnCoefs) V[i] = var(F[,i])
    seqVect = 1:rnCoefs
    indices = seqVect[rank(V) > (rnCoefs-nbCoefs)]

    result = matrix(nrow=length(indices),ncol=ncol(data))
    counter = 1
    for (i in indices) {
        result[counter,] = basis[i,]
        counter = counter + 1
    }
    return (result)
}

#generate Fourier basis on an interval (card == nbCoefs)
genFourier = function(data,nbCoefs,withVar=TRUE) {
    D = ncol(data)
    res = matrix(nrow=floor(D/2)+1,ncol=D)
    res[1,] = rep(1.0/(sqrt(D)),D)
    times = 1:D
    for (i in 1 : D%/%4) {
        res[2*i,] = cos( (2*pi*i*times) / (D-1) )
        res[2*i+1,] = sin( (2*pi*i*times) / (D-1) )
    }

    # orthonormalize the basis :
    for (i in 2:nrow(res)) res[i,] = res[i,] / sqrt(sum(res[i,]^2))
    if (!withVar) return (res[1:min(nbCoefs,nrow(res)),])

    #get the max variance coefficients according to data :
    return (basisMostVar(res,data,nbCoefs))
}

#return the median based on L2 norm ; this choice could be improved..
getMedElem = function(data) {
    n = nrow(data) ; seqVect = 1:n
    nd = sqrt( rowSums(data^2) )
    R = rank(nd)
    med = data[seqVect[R == n%/%2],]
    return (med)
}

# "normal" wavelet basis expansion
simpleWavBasis = function(data, lvl, filt) {
    D = ncol(data)
    med = getMedElem(data)

    W <- wavDWPT(med, n.level = lvl, wavelet = filt)
    wm = wavMRD(W)
    nbFuncs = ncol(wm)
    basis = matrix(nrow=nbFuncs, ncol=D)
    for (i in 1:nbFuncs) basis[i,] = wm[[i]]

    # orthonormalize the basis :
    for (i in 1:nbFuncs) basis[i,] = basis[i,] / sqrt(sum(basis[i,]^2))
    return (basis)
}

#decomposition onto an orthonormal basis :
linEmb = function(data, dim, linbt="PCA",filt="haar",wvar=TRUE) {
    D = ncol(data)

    if (linbt=="PCA") {
        # (functional..) PCA analysis :
        s = svd(data,nu=0,nv=dim)
        basis = t(s$v)
    } else

    if (linbt=="four") {
        if (dim >= ceiling(D/2)) return ( linEmb(data,dim) )
        #Fourier basis :
        basis = genFourier(data,dim,wvar)
    } else

    if (linbt=="wav") {
        lvl = floor(log2(D))
        if (2^lvl < dim) return ( linEmb(data,dim) )
        s = simpleWavBasis(data,lvl,filt)

        basis = matrix()
        if (wvar) basis = basisMostVar(s,data,dim)
        else basis = s[1:min(nrow(s),dim),]
    } else

    if (linbt=="bsp") {
        tmp = t(ns(1:D,df=dim))
        # orthonormalize the basis :
        s = svd(tmp,nu=0,nv=dim)
        basis = t(s$v)
    }

    return (list("embed"=as.matrix(data %*% t(basis)),"basis"=basis))
}


#############################################################
## reconstruction function (easy in those cases..) :

#coefs. in lines, basis functions in lines too
linear_rec = function(coefs, basis) {
    if (!is.matrix(coefs)) coefs = t(as.matrix(coefs))
    return ( coefs %*% basis )
}

