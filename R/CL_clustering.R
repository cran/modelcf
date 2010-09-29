## clustering of inputs or outputs data

#wrapper for hierarchical clustering, Ward linkage :
phclust = function(dissims, K) {
    hct = hclust(as.dist(dissims),method="ward")
    return ( cutree(hct, K) )
}

#main function :
gtclusts = function(method, data, K, d=min(10, ncol(data)), adn=FALSE, knn=0, symm=TRUE, weight=FALSE, sigmo=FALSE, alpha=1.0) {

    if (knn==0) knn = getKnn(nrow(data))

    dissims = matrix()
    if (method=="HDC") dissims = hitorct(data, adn, d, knn, FALSE, symm, weight)
    else if (method=="CTH" || method=="CTHC" || method=="CTKM") dissims = ctdists(data, adn, d, knn, symm, weight, sigmo)
    else if (method=="CH" || method=="CHC") dissims = as.matrix(dist(data))

    clusts = c()
    if (method=="HDC" || method=="CTKM") clusts = km_dists(dissims, K)
    else if (method=="CTHC" || method=="CHC") clusts = chameleon(data, dissims, K, alpha)
    else if (method=="CTH" || method=="CH") clusts = phclust(as.dist(dissims), K)

    else if (method=="specH" || method=="specKM") clusts = spec_clust(method, data, K, d, adn, knn)

    else if (method=="PCA") {
        #always the simplified version (gives good results); the other one is too costly..
        clusts = km_PCA(data, K, d, TRUE)
    }
    else if (method=="KM") clusts = (kmeans(data,K,iter.max=100,nstart=10))$cluster

    return (clusts)
}

