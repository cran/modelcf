############################################
##selective libraries loading :
checkDependP = function(mclust="",mclass="",mdim="",linbt="",mreg="") {
    if (mclust=="spec") require("kernlab",quietly=TRUE)
    if (mclass=="kNN") require(class,quietly=TRUE)
    if (mclass=="SVM" || mreg=="SVR") require(e1071,quietly=TRUE)
    if (mclass=="ctree") require(rpart,quietly=TRUE) #classification tree
    if (mclass=="rforest" || mreg=="rforest") require(randomForest,quietly=TRUE)
    if (mclass=="RDA") require(klaR,quietly=TRUE) #Regularized Discriminant Analysis
    if (mreg=="GP") require(mlegp,quietly=TRUE)
    if (linbt=="wav") require(wmtsa,quietly=TRUE) #wavelets
    if (linbt=="bsp") require(splines,quietly=TRUE) #B-splines
}


############################################
## useful for testing/building metamodel :

#return a curve matching the maximas given by user (to facilitate model mixing)
getcoefc = function(D, inds, maxs=c(), rgs=c()) {
    res = rep(0.0,D)
    nbInds = length(inds)
    if (length(maxs) == 0) {
        #set all maxs to 1 :
        maxs = rep(1, nbInds)
    }
    if (length(rgs) == 0) {
        #look for the smallest gap
        smallestGap = inds[1]-1
        for (i in 2:nbInds) {
            delta = inds[i] - inds[i-1]
            if (delta < smallestGap) smallestGap = delta
        }
        rgs = rep(smallestGap, nbInds)
    }
    for (i in 1:nbInds) {
        delta = round(rgs[i]/2.0)
        res[ max(1, inds[i]-delta) : min(D, inds[i]+delta) ] = maxs[i]
    }
    return (res)
}

#compute MSE and Q2 errors indicators, by comparing predictions to real curves :
fperrors = function(ypred, yreal, mntrain) {

    #if arguments are file names :
    if (is.character(ypred)) ypred = as.matrix(read.table(ypred))
    if (is.character(yreal)) yreal = as.matrix(read.table(yreal))

    #if only one output is given :
    if (!is.matrix(ypred)) ypred = t(as.matrix(ypred))
    if (!is.matrix(yreal)) yreal = t(as.matrix(yreal))

    nbEx = nrow(ypred) ; D = ncol(ypred)
    MSE = matrix(0.0,nrow=nbEx,ncol=D)
    pvar = matrix(0.0,nrow=nbEx,ncol=D)
    for (i in 1:nbEx) {
         MSE[i,] = (yreal[i,] - ypred[i,])^2
         pvar[i,] = (yreal[i,] - mntrain)^2
    }
    return ( list("MSE"=colMeans(MSE),"pvar"=colMeans(pvar)) )
}


############################################
## MISC :

#standardize data :
standz = function(data) {
    mean = colMeans(data)
    stdevs = sqrt( unlist( apply(data, 2, var) ) )
    res = t(data) - mean
    res = t(res / stdevs)
    return (list("data"=res,"mean"=mean,"stdevs"=stdevs))
}

#Moore-Penrose pseudo inverse
mppsinv = function(M) {
    s = svd(M)
    sd = s$d ; sd[sd < EPS()] = Inf
    sd = diag(1.0 / sd, min(nrow(M),ncol(M)))
    return ((s$v) %*% sd %*% t(s$u))
}

#if simple knn needed : empirical heuristic..
getKnn = function(n) {
    return ( max(1, min( n-1, ceiling(sqrt(n)) ) ) ) #+ floor(log(n)/log(5)) ) ) )
}


############################################
## PLOTTING and PRINTING :

#plot a matrix of curves (in rows)
plotC = function(data, cl=rep(1,nrow(data)), rg=c(min(data),max(data)), ...) {
    n = nrow(data)
    for (i in 1:n) {
        plot(data[i,],col=cl[i],ylim=rg,type="l", ...)
        if (i < n) par(new=TRUE)
    }
}

#for scatterplots : easier than lattice..
plotPts = function(data, cols=c(1,2), cl=rep(1,nrow(data)), ...) {
    rg1=range(data[,cols[1]]) ; rg2=range(data[,cols[2]])
    n = nrow(data)
    for (i in 1:n) {
        plot(data[i,cols[1]],data[i,cols[2]],pch=16,col=cl[i],xlim=rg1,ylim=rg2, ...)
        if (i < n) par(new=TRUE)
    }
}

#simple print method
print.modelcf = function(x, ...) {
    model = x
    n = nrow(model$x) ; p = ncol(model$x)
    D = ncol(model$y) ; d = model$d
    nbCs = max(model$clusts)
    cat(paste("\nclass of object: ",pkgnm(),"\n",sep=""))
    cat(paste("model built on ",n," samples\n",sep=""))
    cat(paste("    inputs / outputs dimensionality: ",p," / ",D,"\n",sep=""))

    if (nbCs>1) {
        cat("\n=========================\n")
        cat(paste("\n",nbCs," clusters found; use plotC(y,cl=model$clusts) to visualize\n",sep=""))
        type = model$mclust ; clname = ""
        if (type=="HDC") clname="k-means based on Hitting Times"
        else if (type=="CTH") clname="Commute-Time Hierarchic"
        else if (type=="CTKM") clname="Commute-Time k-means"
        else if (type=="specH") clname="spectral clustering (with hierarchical clust.)"
        else if (type=="specKM") clname="spectral clustering (with k-means)"
        else if (type=="CH") clname="hierarchical clustering"
        else if (type=="ACP") clname="ACP-k-means from Chiou and Li"
        else if (type=="KM") clname="basic k-means"
        cat(paste("    inputs-outputs clustering method: ",clname,"\n",sep=""))
        type = model$mclass ; csname = ""
        if (type=="kNN") csname="k Nearest Neighbors"
        else if (type=="ctree") csname="classification trees"
        else if (type=="SVM") csname="Support Vector Machines"
        else if (type=="rforest") csname="Random Forests"
        else if (type=="RDA") csname="Regularized Discriminant Analysis"
        cat(paste("    inputs classification method: ",csname,"\n",sep=""))
    }

    if (d<D) {
        cat("\n=========================\n")
        cat(paste("\noutputs reduced to ",d," dimensions, using\n",sep=""))
        type = model$mdim ; dname = ""
        if (type=="linear") {
            dname = "orthogonal basis"
            linbt = model$linbt
            if (linbt == "PCA") dname = paste(dname,": functional PCA",sep="")
            else if (linbt == "four") dname = paste(dname,": Fourier basis",sep="")
            else if (linbt == "bsp") dname = paste(dname,": B-splines basis",sep="")
            else if (linbt == "wav") dname = paste(dname,": wavelets basis, filter= ",model$filt,sep="")
        }
        else if (type=="RML") dname="Riemannian Manifold Learning"
        else if (type=="LPcaML") dname="Local PCA Manifold Learning"
        else if (type=="GCEM") dname="Local RML Manifold Learning"
        cat(paste("    ",dname,"\n",sep=""))
    }

    cat("\n=========================\n")
    type = model$mreg ; rname = ""
    if (type=="PPR") rname="Projection Pursuit Regression"
    else if (type=="BRT") rname="Boosting Regression Trees"
    else if (type=="SVR") rname="Support Vector Regression"
    else if (type=="rforest") rname="Random Forests"
    else if (type=="GP") rname="Gaussian Processes"
    else if (type=="kNN") rname="Nadaraya-Watson on reduced representations"
    else if (type=="fkNN") rname="Nadaraya-Watson directly on outputs"
    cat(paste("\nregression method: ",rname,"\n",sep=""))
    if (model$stred) cat("    standardized outputs\n")
    if (model$ppts) cat("    pointwise regression\n")

    cat("\n")
}

