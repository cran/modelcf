#compare two partitions...

## (simple) assymmetric measure :

#intuition : for each matching tentative, we compute
#delta = max. clusters' overlap - second max. cluster's overlap.
checkParts = function(P, P_ref) {
    n = length(P_ref) ; seqVect = 1:n
    accuracy = 0
    K = max(P, P_ref)
    for (i in 1:K) {
        #we will compare cluster i of P..
        clustI_P = (P == i)
        firstMax = 0 ; secMax = 0
        for (j in 1:K) {
            #..to cluster j of P_ref :
            clustJ_Pref = (P_ref == j)
            tmp = sum( clustI_P & clustJ_Pref ) / sum(clustJ_Pref)

            if (tmp > firstMax) {
                secMax = firstMax
                firstMax = tmp
            }
            else if (tmp > secMax) secMax = tmp
        }
        accuracy = accuracy + firstMax - secMax
    }
    return (accuracy / K)
}


## symmetric measures :

#Meila 2002 :
varInfo = function(P1, P2) {
    n = length(P1) ; seqVect = 1:n
    K = max(P1, P2)

    #some pre-calculations :
    clusts1 = as.list(rep(0,K))
    clusts2 = as.list(rep(0,K))
    Pk1 = rep(0.0,K)
    Pk2 = rep(0.0,K)
    for (i in 1:K) {
        clusts1[[i]] = seqVect[P1 == i]
        clusts2[[i]] = seqVect[P2 == i]
        Pk1[i] = length(clusts1[[i]]) / n
        Pk2[i] = length(clusts2[[i]]) / n
    }

    #a few lines to get VI :
    Pk1_ = Pk1[Pk1 > 0.0]
    Pk2_ = Pk2[Pk2 > 0.0]
    entrop1 = 0.0
    entrop2 = 0.0
    if (length(Pk1_) > 0) entrop1 = - sum( Pk1_ * log(Pk1_) )
    if (length(Pk2_) > 0) entrop2 = - sum( Pk2_ * log(Pk2_) )
    I = 0.0
    for (i in 1:K) {
        for (j in 1:K) {
            Pkk = length( intersect(clusts1[[i]],clusts2[[j]]) ) / n
            if ( Pkk > 0.0 )
                I = I + Pkk * log( Pkk / (Pk1[i] * Pk2[j]) )
        }
    }
    VI = entrop1 + entrop2 - 2*I
    return (1.0 - VI / log(n) )
}

#simple count index :
countPart = function(P1,P2) {
    n = length(P1)
    res = 0
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            if ( (P1[i]==P1[j] && P2[i]==P2[j])
               || (P1[i]!=P1[j] && P2[i]!=P2[j]) ) res = res + 1
        }
    }
    return ( 2.0 * res / (n * (n-1) ) )
}

