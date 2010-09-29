#######################################
## Learning and prediction :

# method to learn a regression relationship y~x
learnRegress = function(x, y, method, params, stred) {
    x = as.matrix(x) ; y = as.matrix(y)
    d = ncol(y)

    #projection pursuit regression
    if (method=="PPR") {
        sdins = standz(x)
        outs = y ; sdouts = NULL
        if (stred) {
            sdouts = standz(y)
            outs = sdouts$data
        }

        pp = as.list(rep(0,d))
        for (i in 1:d) pp[[i]] = ppr( sdins$data,outs[,i],nterms=params[i] )
        return ( list("type"=method,"obj"=pp,"sdins"=sdins,"sdouts"=sdouts) )
    } else

    #random forest
    if (method=="rforest") {
        outs = y ; sdouts = NULL
        if (stred) {
            sdouts = standz(y)
            outs = sdouts$data
        }

        ra = as.list(rep(0,d))
        for (i in 1:d) ra[[i]] = randomForest(x, outs[,i])
        return ( list("type"=method,"obj"=ra,"sdouts"=sdouts) )
    } else

    # boosting of regression trees :
    if (method=="BRT") {
        outs = y ; sdouts = NULL
        if (stred) {
            sdouts = standz(y)
            outs = sdouts$data
        }

        dtf = data.frame(x)
        tr = as.list(rep(0,d))
        for (i in 1:d) tr[[i]] = gbm(outs[,i] ~ .,data=dtf,distribution="gaussian",n.trees=params[i])
        return ( list("type"=method,"btree"=tr,"ntrees"=params,"df"=dtf,"sdouts"=sdouts) )
    } else

    # Nadaraya-Watson (with or without dim. red.), or local PCA regression
    if (method=="kNN" || method=="fkNN" || method=="lPCA") {
        sdins = standz(x)
        return ( list("type"=method,"x_train"=sdins$data,"y_train"=y,"knn"=params,"sdins"=sdins) )
    } else

    # gaussian processes (quite slow)
    if (method == "GP") {
        sdins = standz(x)
        outs = y ; sdouts = NULL
        if (stred) {
            sdouts = standz(y)
            outs = sdouts$data
        }

        mlg = as.list(rep(0,d))
        for (i in 1:d) mlg[[i]] = mlegp(sdins$data,outs[,i],constantMean=0,seed=sample(1:1000,1))
        return ( list("type"=method,"obj"=mlg,"sdins"=sdins,"sdouts"=sdouts) )
    } else

    # support vector machine (parameters optimisation is very slow)
    if (method=="SVR") {
        outs = y ; sdouts = NULL
        if (stred) {
            sdouts = standz(y)
            outs = sdouts$data
        }

        sv = as.list(rep(0,d))
        for (i in 1:d)
            sv[[i]] = svm(x, outs[,i], type = "eps-regression", kernel = "radial",
                          gamma = params[1,i], cost = params[2,i], epsilon = params[3,i])
        return (list("type"=method,"obj"=sv,"sdouts"=sdouts) )
    }

    return (NULL)
}

#prediction :
predictRegress = function(model,newIn_s) {

    if (!is.matrix(newIn_s)) newIn_s = t(as.matrix(newIn_s))
    #scale input(s) if needed :
    if (length(model$sdins) > 0) newIn_s = t( ( t(newIn_s) - model$sdins$mean ) / model$sdins$stdevs )
    nb2pred = nrow(newIn_s) ; pr = matrix(nrow=nb2pred,ncol=0)

    if (model$type=="BRT") {
        dtf = model$df
        dtf[1:nb2pred,] = newIn_s
        d = length(model$btree)
        for (i in 1:d) pr = cbind(pr, as.matrix(predict(model$btree[[i]],dtf[1:nb2pred,],n.trees=model$ntrees[i])))
    } else

    if (model$type=="kNN" || model$type=="fkNN") { #Nadaraya-Watson (on reduced coordinates or functions)
        pr = knnPredict(model$x_train, model$y_train, newIn_s, model$knn)
    } else

    if (model$type=="lPCA") { #local PCA regression (on functions) :
        pr = lpcaPredict(model$x_train, model$y_train, newIn_s, model$knn)
    }

    else { #all other cases
        d = length(model$obj)
        for (i in 1:d) pr = cbind(pr, as.matrix( predict(model$obj[[i]], newIn_s) ))
    }

    #descale outputs if needed :
    if (length(model$sdouts) > 0) pr = t( t(pr) * (model$sdouts$stdevs) + model$sdouts$mean )
    if (nrow(pr)==1) pr = as.double(pr)
    return (pr)
}


#######################################
## Parameters optimisation :

#parameters optimisation with a "good" test set :
optimParams_regress = function(x,y,method,knn,trcv,verb) {
    if (method=="rforest" || method=="GP") return (0)

    #extract training data from x rows :
    design = xtr_plan2(x,knn,trcv)
    x_train = as.matrix(x[design,])
    x_test = as.matrix(x[-design,])
    y_train = as.matrix(y[design,])
    y_test = as.matrix(y[-design,])

    #and optimize the chosen model parameters on those sets : [N-W and lPCA = special]
    if (method=="kNN" || method=="fkNN" || method=="lPCA") {
        bestInd=1 ; bestObj = Inf
        for (j in 1:( min( nrow(x_train)%/%5, 50 ) )) {
            kl_obj = learnRegress(x_train,y_train,method,j,FALSE)
            pr = predictRegress(kl_obj, x_test)
            op = sum( sqrt( rowSums( (pr - y_test)^2 ) ) )
            if (op < bestObj) {
                bestInd = j
                bestObj = op
            }
        }
        return (bestInd)
    }

    res = 0 ; d = ncol(y)
    if (method=="SVR") res = matrix(nrow=3,ncol=d)
    else if (method=="PPR") res = rep(0,d)

    for (i in 1:d) {
        if (verb) cat(paste("    optimize parameters for dimension ",i,"\n",sep=""))

        if (method=="SVR") {
            tc = tune.control()
            tc$cross = min(nrow(x_train),10)
            svrParams = tune(svm, train.x=x_train, train.y=y_train[,i], 
                validation.x = x_test, validation.y = y_test[,i],
                ranges = list(gamma = 2^(-2 : 1),cost = 2^(-1 : 2), epsilon = 2^(-4 : -1)),
                tunecontrol = tc, type="eps-regression")
            res[1,i] = svrParams$best.parameters$gamma
            res[2,i] = svrParams$best.parameters$cost
            res[3,i] = svrParams$best.parameters$epsilon
        } else

        if (method=="PPR" || method=="BRT") {
            bestInd=1 ; bestObj = Inf
            for (j in 1:(min(3*d,30))) {
                learn_obj = learnRegress(x_train,y_train[,i],method,j,TRUE)
                pr = predictRegress(learn_obj,x_test)
                op = sqrt( sum( (pr - y_test[,i])^2 ) )
                if (op < bestObj) {
                    bestInd = j
                    bestObj = op
                }
            }
            res[i] = bestInd
        }
    }

    return (res)
}

