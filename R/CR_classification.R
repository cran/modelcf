#######################################
## Learning and prediction :

# method to build a classification object
learnClassif = function(x, y, method, params) {
    x = as.matrix(x)

    #k nearest neighbors
    if (method=="kNN") {
        sdins = standz(x) #use standardized data
        return ( list("type"=method,"x_train"=sdins$data,"y_train"=y,
                      "meanx"=sdins$mean,"stdevsx"=sdins$stdevs,"knn"=params) )
    } else

    #classification tree
    if (method=="ctree") {
        df = data.frame(x)
        tem = df[1,]
        tr = rpart(y ~ ., data=df, method="class")
        return ( list("type"=method,"tree"=tr,"tm"=tem) )
    } else

    # Regularized Discriminant Analysis
    if (method=="RDA") {
        sdins = standz(x) #use standardized data
        res = rda(sdins$data, grouping=y, train.fraction=0.7, estimate.error=FALSE)
        return ( list("type"=method,"obj"=res,"meanx"=sdins$mean,
                      "stdevsx"=sdins$stdevs) )
    } else

    #random forest
    if (method=="rforest") {
        res = randomForest(x, as.factor(y))
        return ( list("type"=method,"obj"=res) )
    } else

    # support vector machines
    if (method == "SVM") {
        res = svm(x, y, type = "C-classification", kernel = "radial",
                    gamma = params[1], cost = params[2], epsilon = params[3])
        return ( list("type"=method,"obj"=res) )
    }

    return (NULL)
}

#predict method associated :
predictClassif = function(model,newIns) {
    if (!is.matrix(newIns)) newIns = t(as.matrix(newIns))
    clElems = c()

    if (model$type == "kNN") {
        scIns = t( ( t(newIns) - model$meanx ) / model$stdevsx )
        clElems = as.integer( knn(train=model$x_train, test=scIns, cl=model$y_train, k=model$knn) )
    }

    else if (model$type=="ctree") {
        tem = model$tm
        for (i in 1:nrow(newIns)) {
            tem[] = newIns[i,]
            p = predict(model$tree, tem)
            clElems = c(clElems,which.max(p))
        }
    }

    else if (model$type=="RDA") {
        scIns = t( ( t(newIns) - model$meanx ) / model$stdevsx )
        clElems = predict( model$obj, scIns )
        clElems = as.integer(p[[1]])
    }

    else clElems = predict( model$obj, newIns )
    return (clElems)
}


#######################################
## Parameters optimisation :

#parameters optimization with a "good" test set :
optimParams_classif = function(x,y,method,knn,trcv) {
    if (method!="kNN" && method!="SVM") return (0)

    #extract training data from x rows :
    design = xtr_plan2(x,knn,trcv)
    x_train = as.matrix(x[design,])
    x_test = as.matrix(x[-design,])
    y_train = y[design]
    y_test = y[-design]

    #optimize the chosen model parameters on these sets :
    res = 0
    if (method=="SVM") {
        svmParams = tune(svm, train.x=x_train, train.y=y_train, ranges = list(gamma = 2^(-2:1),
                         cost = 2^(-1 : 2), epsilon = 2^(-4 : -1)), type="C-classification")
        res = rep(0.0,3)
        res[1] = svmParams$best.parameters$gamma
        res[2] = svmParams$best.parameters$cost
        res[3] = svmParams$best.parameters$epsilon
    }
    else if (method=="kNN") {
        bestInd=1 ; bestObj = Inf
        for (j in 1:max(1, min(30, nrow(x_train)%/%10) ) ) {
            cl = as.integer( knn(train=x_train, test=x_test, cl=y_train, k=knn) )
            op = sum( (cl - y_test) & 1)
            if (op < bestObj) {
                bestInd = j
                bestObj = op
            }
        }
        res = bestInd
    }
    return (res)
}

