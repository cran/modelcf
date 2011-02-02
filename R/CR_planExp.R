#extract a training set and test set from n samples
#..using local variance :
xtr_plan1 = function(data, k, trcv) {

    #reduce dimension if too many :
    if (ncol(data) > 4) {
        s = svd(data,nu=0,nv=4)
        data = (data %*% s$v) %*% t(s$v)
    }
    data = (standz(data))$data #..and standardize

    #compute all local variances (inputs ==> no outliers) :
    n = nrow(data) ; locVar = rep(0.0,n)
    for (i in 1:n) {
        dists = sort( rowSums( ( t( t(data) - data[i,]) )^2 ) )
        locVar[i] = mean(dists[2:(k+1)])
    }
    srt = sort(locVar,index.return=TRUE)

    # finally get the design from these local variances :
    nbtrain = min( n-2, max(2,floor(trcv*n)) )
    return (srt$ix[1:nbtrain])
}

#..using local densities :
xtr_plan2 = function(data, k, trcv) {

    #reduce dimension if too many :
    if (ncol(data) > 4) {
        s = svd(data,nu=0,nv=4)
        data = data %*% s$v
    }
    data = (standz(data))$data #..and standardize

    #compute all distances and sort them to obtain dist. to knn neighb.
    n = nrow(data) ; locDists = rep(0.0,n)
	for (i in 1:n) {
        dists = sort( rowSums( ( t( t(data) - data[i,]) )^2 ) )
        locDists[i] = dists[k+1]
    }
	srt = sort(locDists, index.return=TRUE)

    # finally get the design from these local densities :
    nbtrain = min( n-2, max(2,floor(trcv*n)) )
    return (srt$ix[1:nbtrain])
}

