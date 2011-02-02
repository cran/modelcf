#######################################
## Learning and prediction :

# method to build a classification object
learnClassif = function(x, y, method, params) {
    x = as.matrix(x)
    checkDependP(mclass=method)

    #k nearest neighbors
    if (method=="kNN") {
        sdins = standz(x) #use standardized data
        return ( list("type"=method,"x_train"=sdins$data,"y_train"=y,
                      "meanx"=sdins$mean,"stdevsx"=sdins$stdevs,"kNN"=params) )
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
        clElems = as.integer( knn(train=model$x_train, test=scIns, cl=model$y_train, k=model$kNN) )
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
optimParams_classif = function(x,y,method,k,trcv) {
	checkDependP(mclass=method)
    if (method!="kNN" && method!="SVM") return (0)

	design = c()
	n = nrow(x)
	if (n < 10) {
		design = 1:n
    	x_test = as.matrix(x[design,])
		y_test = as.matrix(y[design])
    }
    #extract training data from x rows :
    else {
		design = xtr_plan2(x,k,trcv)
		x_test = as.matrix(x[-design,])
		y_test = as.matrix(y[-design])
    }
    x_train = as.matrix(x[design,])
    y_train = as.matrix(y[design])

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
        bestK=1 ; bestObj = Inf
        nTR = nrow(x_train)
        knnTR = getKnn(nTR)
        for (j in 1:min(30, 2*knnTR, nTR-1) ) {
			cl = as.integer( knn(train=x_train, test=x_test, cl=y_train, k=j) )
			op = sum( (cl - y_test) & 1)
			if (op < bestObj) {
				bestK = j
				bestObj = op
			}
		}
		res = bestK
    }
    return (res)
}

