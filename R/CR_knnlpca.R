#Nadaraya-Watson estimate :
knnPredict = function(x_train, y_train, x, k) {
    m = nrow(x)
    pr = matrix(0.0,nrow=m,ncol=ncol(y_train))

    for (i in 1:m) {
        cs = colSums( ( t(x_train) - x[i,])^2 )
        srt = sort(cs, index.return=TRUE)
        if ( k==1 || sqrt(srt$x[1]) <= EPS() ) pr[i,] = y_train[srt$ix[1],]
        else {
            sigma2 = (srt$x[k] - srt$x[1]) / log(srt$x[k] / srt$x[1])
            coefs = exp(- srt$x[1:k] / sigma2)
            coefs = coefs / sum(coefs)
            if (k > 1) pr[i,] = coefs %*% y_train[srt$ix[1:k],]
            else pr[i,] = coefs * y_train[srt$ix[1],]
        }
    }
    return (pr)
}

# kind of local PCA regression
lpcaPredict = function(x_train, y_train, x, kd, stzouts=TRUE) {
    m = nrow(x)
    n = nrow(y_train)
    pr = matrix(0.0,nrow=m,ncol=ncol(y_train))
    k = kd[1]
    d = kd[2]

    for (i in 1:m) {
        #1) PCA around data[nearest_neighbor_of_newEmb,] :
        cs = colSums( (t(x_train) - x[i,])^2 )
        srt = sort(cs, index.return=TRUE)
        if ( k==1 || sqrt(srt$x[1]) <= EPS() ) pr[i,] = y_train[srt$ix[1],]
        else {
        
			nbEx = min(n-1, 2*k) #number of examples for machine learning
			inds = srt$ix[1:(nbEx+1)]

			#get the LOCAL basis, and the coefficients in the neighborhood :
			inds_SVD = srt$ix[1:(max(d+1,k))]
			s = svd(y_train[inds_SVD,],nu=0,nv=d)
			basis = as.matrix(s$v)
			allACP_coefs = y_train[inds,] %*% basis

			#2) estimate optimal learning parameters
			#no time here for slow optimization ; so, fast heuristic :
			op = rep( min(d, 30), d )

			#3) "solve" by statistical learning on the coefficients VS. PCA coefs.
			regCoefs = learnRegress(x_train[inds,], allACP_coefs, "PPR", op, stzouts)
			res = predictRegress(regCoefs, x[i,])

			#4) linear reconstruction
			pr[i,] = res %*% t(basis)
		}
    }

    return (pr)
}

