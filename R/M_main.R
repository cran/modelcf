#######################################
##training methods for the metamodel :

#dimension reduction y_i --> z_i and regression x_i --> z_i steps (within one cluster)
fmetam_1cl = function(x, y, d,
                    ## params red dim :
                    mdim,adnRD,knnRD,linbt,filt,wvar,alpha,knnmin,knnmax,tsoft,thlvl,hdth,
                    ## params regression :
                    mreg, ppts, stred, trcv, verb) {

    n = nrow(y)
    embobj = NULL
    if (!ppts && mreg!="fkNN" && mreg!="lPCA") {
        if (verb) cat("\n*** dimension reduction ***\n")
        if (mdim == "linear") embobj = linEmb(y, d, linbt, filt, wvar)
        else embobj = nlin_redDim(mdim,y,d,adnRD,knnRD,alpha,knnmin,knnmax,tsoft,thlvl,hdth)
    }
    else embobj = list("embed"=y)

    if (verb) cat("\noptimize regression parameters\n")
    regressParams = optimParams_regress(x, embobj$embed, mreg, getKnn(n), trcv, verb)
    if (verb) cat("\n*** regression ***\n")
    outModel = learnRegress(x, embobj$embed, mreg, regressParams, stred)

    return (list("model"=outModel,"embobj"=embobj))
}

#main call : metamodeling
fmetam = function(x, y, d=0,
                    ## params clust et classif :
                    mclust="CTH",mclass="kNN",redy=TRUE,adnCC=TRUE,knnCC=0,wcl=TRUE,symm=FALSE,
                    weight=TRUE,sigmo=FALSE,alphaCL=1.0,minszcl=30,maxcl=Inf,taus=0.95,Ns=10,tauc=0.95,Nc=10,
                    ## params red dim :
                    mdim="linear",adnRD=FALSE,knnRD=0,linbt="PCA",filt="haar",wvar=TRUE,alphaRD=0.5,
                    knnmin=0,knnmax=0,tsoft=0.1,thlvl=0.3,hdth=0,
                    ## params regression :
                    mreg="PPR", ppts=FALSE, stred=TRUE, trcv = 0.7, verb = TRUE) {

    #if arguments are file names :
    if (is.character(x)) x = as.matrix(read.table(x))
    if (is.character(y)) y = as.matrix(read.table(y))

    n = nrow(y) ; D = ncol(y)
    p = ncol(x) ; seqVect = 1:n
    #estimate dimension if needed :
    if (ppts) d = D
    else if (mreg!="fkNN" && d==0) {
        d = min( p, dimest1(y, getKnn(n) ) )
        if (verb) cat(paste("\nestimated dimension = ",d,"\n",sep=""))
    }

    #check dependencies / load appropriate libraries :
    checkDependP(mreg,mclass,mdim,linbt)

    #####################################################

    clusts = rep(1,n)
    if (wcl) {
        if (verb) cat("\n*** clustering ***\n")
        clusts = gtclusts_inout(x,y,mclust,d,redy,adnCC,knnCC,symm,weight,sigmo,alphaCL,minszcl,maxcl,mclass,taus,Ns,tauc,Nc,trcv,verb)
    }
    nbClusts = max(clusts)
    if (verb && wcl) {
        if (nbClusts > 1) cat(paste(nbClusts," clusters have been found\n",sep=""))
        else cat("no clusters have been found\n")
    }

    classifParams = 0
    modclass = NULL
    if (nbClusts > 1) {
        if (verb) cat("\noptimize classification parameters\n")
        classifParams = optimParams_classif(x,clusts,mclass,knn,trcv)
        if (verb) cat("\n*** classification ***\n")
        modclass = learnClassif(x,clusts,mclass,classifParams)
    }

    #####################################################

    modsreg = as.list( rep(0,nbClusts) )
    for (i in 1:nbClusts) {
        if (verb) cat(paste("\n=== cluster ",i," with ",sum(clusts==i)," elements: ===\n"))
        modsreg[[i]] = fmetam_1cl(as.matrix(x[clusts == i,]),y[clusts == i,],d,mdim,adnRD,knnRD,linbt,
                           filt,wvar,alphaRD,knnmin,knnmax,tsoft,thlvl,hdth,mreg,ppts,stred,trcv,verb)
    }

    if (verb) cat("\n")
    return ( structure( list("x"=x,"y"=y,"d"=d,"mclust"=mclust,"mclass"=mclass,"redy"=redy,"adnCC"=adnCC,"knnC"=knnCC,
                             "clusts"=clusts,"symm"=symm,"weight"=weight,"sigmo"=sigmo,"alphaCL"=alphaCL,"minszcl"=minszcl,"maxcl"=maxcl,
                             "modclass"=modclass,"taus"=taus,"Ns"=Ns,"tauc"=tauc,"Nc"=Nc,"mdim"=mdim,"adnRD"=adnRD,
                             "knnRD"=knnRD,"linbt"=linbt,"filt"=filt,"wvar"=wvar,"alphaRD"=alphaRD,"knnmin"=knnmin,
                             "knnmax"=knnmax,"tsoft"=tsoft,"thlvl"=thlvl,"hdth"=hdth,"mreg"=mreg,"modsreg"=modsreg,
                             "ppts"=ppts,"stred"=stred,"trcv"=trcv,"verb"=verb), class=pkgnm() ) )
}


#######################################
## testing methods for the metamodel :

#..with a new test matrix (or vector) :
predict.modelcf = function(object, x, verb = FALSE, ...) {

    #if argument is file name :
    if (is.character(x)) x = as.matrix(read.table(x))
    #if only one output is given
    if (!is.matrix(x)) x = t(as.matrix(x))

    mo = object ; d = mo$d
    y_train = mo$y ; D = ncol(y_train)
    nbEx = nrow(x)
    predCurves = matrix(nrow=nbEx,ncol=D)

    #stage 1: find classes of testing elements, following mo$inClassif :
    clElems = rep(1,nbEx)
    if ( max(mo$clusts) > 1 ) {
        if (verb) cat("\n*** classify testing data ***\n")
        clElems = predictClassif(mo$modclass,x)
    }

    #stage 2: predict reduced coordinates for each test vector
    for (i in 1:nbEx) {

        if (verb) {
            cat(paste("\ntest ",i,"th input\n",sep=""))
            cat("    *** predict (reduced) coordinates ***\n")
        }
        pr = predictRegress( (mo$modsreg[[ clElems[i] ]])$model, t(as.matrix(x[i,])) )
        rec = rep(0.0,D)

        if (mo$mreg == "fkNN" || mo$mreg == "lPCA" || mo$ppts)
            #Nadaraya-Watson, local PCA or pointwise :
            rec = pr

        else {
            #other cases :
            if (verb) cat("    *** build curve from reduced representation ***\n")

            if (mo$mdim == "linear") rec = linear_rec(pr, (mo$modsreg[[ clElems[i] ]])$embobj$basis )

            else {
                #method adaptive (for RML, LPcaML and GCEM -- others are UNAVAILABLE for metamodeling,
                #because they haven't adapted reconstruction algorithms)
                rec = nlin_adaptRec(mo$mdim, mo$modsreg[[ clElems[i] ]]$embobj, pr)
            }
        }

        predCurves[i,] = rec
    }

    if (verb) cat("\n")
    return (predCurves)
}

#prediction for a mixture of models :
mixpredf = function(mods, coefs, x, verb=FALSE) {

    #if argument is file name :
    if (is.character(x)) x = as.matrix(read.table(x))
    #if only one output is given
    if (!is.matrix(x)) x = t(as.matrix(x))

    nbMods = length(mods)
    D = length(coefs[[1]])
    #normalize coefficient curves to sum to 1 at each point when possible :
    for (i in 1:D) {
        psum = 0.0
        for (j in 1:nbMods) psum = psum + coefs[[j]][i]
        if ( psum > EPS() ) {
            for (j in 1:nbMods) coefs[[j]][i] = coefs[[j]][i] / psum
        }
    }
    preds = matrix(0.0, nrow=nrow(x), ncol=D)
    for (i in 1:nbMods) {
		locPreds = predict.modelcf(mods[[i]], x, verb)
		preds = preds + t(coefs[[i]] * t(locPreds))
	}
    return (preds)
}

#generalized cross-validation
nfoldcv = function(x, y, d=0,
                    ## params clust et classif :
                    mclust="CTH",mclass="kNN",redy=TRUE,adnCC=TRUE,knnCC=0,wcl=TRUE,symm=FALSE,
                    weight=TRUE,sigmo=FALSE,alphaCL=1.0,minszcl=30,maxcl=Inf,taus=0.95,Ns=10,tauc=0.95,Nc=10,
                    ## params red dim :
                    mdim="linear",adnRD=FALSE,knnRD=0,linbt="PCA",filt="haar",wvar=TRUE,alphaRD=0.5,
                    knnmin=0,knnmax=0,tsoft=0.1,thlvl=0.3,hdth=0,
                    ## params regression :
                    mreg="PPR", ppts=FALSE, stred=TRUE, trcv = 0.7,
                    ## params cross-validation and printing options :
                    loo = FALSE, nfold=100, nhold=10, verb = TRUE, plotc=TRUE) {

    #if arguments are file names :
    if (is.character(x)) x = as.matrix(read.table(x))
    if (is.character(y)) y = as.matrix(read.table(y))

    n = nrow(y) ; D = ncol(y) ; seqVect = 1:n
    MSE = matrix(nrow=nfold,ncol=D)
    pvar = matrix(nrow=nfold,ncol=D)
    Q2 = matrix(nrow=nfold,ncol=D)
    statsNbClusts = rep(0,n) ; sameSizeClusts = 0

    predCurves = matrix()
    if (loo) {
        nfold = n
        predCurves = matrix(nrow=n,ncol=D)
    }

    for (i in 1:nfold) {
        trainInds = c() ; testInds = c()
        if (loo) {
            trainInds = seqVect[-i]
            testInds = i
        }
        else {
            testInds = sample(seqVect,nhold)
            trainInds = seqVect[-testInds]
        }

        if (verb) cat(paste("\n||| build metamodel ",i," |||\n",sep=""))
        met = fmetam(as.matrix(x[trainInds,]),y[trainInds,],d,mclust,mclass,redy,adnCC,knnCC,wcl,symm,
                     weight,sigmo,alphaCL,minszcl,maxcl,taus,Ns,tauc,Nc,mdim,adnRD,knnRD,linbt,filt,wvar,alphaRD,
                     knnmin,knnmax,tsoft,thlvl,hdth, mreg, ppts, stred, trcv, verb)
        nbClusts = max(met$clusts) ; minSize = Inf
        for (j in 1:nbClusts) {
            tmp = sum(met$clusts==j)
            if (tmp < minSize) minSize = tmp
        }
        sameSizeClusts = sameSizeClusts + max(0, n%/%nbClusts - minSize)
        statsNbClusts[nbClusts] = statsNbClusts[nbClusts]+1

        if (verb) cat(paste("\n||| test metamodel ",i," |||\n",sep=""))
        pcurves = predict.modelcf(met, x[testInds,], verb)
        errs = fperrors(pcurves, y[testInds,], colMeans(y[trainInds,]))
        MSE[i,] = errs$MSE
        pvar[i,] = errs$pvar
        Q2[i,] = 1.0 - MSE[i,]/pvar[i,]
        if (loo) predCurves[i,] = pcurves[1,]

        if (plotc) {
            plot(Q2[i,],ylim=c(0,1),type="l",ylab="prediction Q2",cex.lab=2,cex.axis=2)
            par(new=TRUE)
        }
        if (verb) cat("=========================\n\n")
    }

    #empirical standard deviations on error estimators :
    stdMSE = rep(0.0,D) ; stdPvar = rep(0.0,D) ; stdQ2 = rep(0.0,D)
    for (i in 1:D) {
        stdMSE[i] = sqrt(var(MSE[,i]))
        stdPvar[i] = sqrt(var(pvar[,i]))
        stdQ2[i] = sqrt(var(Q2[,i]))
    }
    MSE = colMeans(MSE)
    pvar = colMeans(pvar)
    Q2 = 1.0 - MSE / pvar

    return ( list("curves" = predCurves, "MSE"=MSE, "stMSE"=stdMSE, "pvar"=pvar, "stvar"=stdPvar,
                  "Q2"=Q2, "stQ2"=stdQ2, "ssclust"=sameSizeClusts/nfold, "snclust"=statsNbClusts) )
}

