###########################################################
## methods to reorganize clustering :

# renumbering step :
reordering = function(clusts) {
    newCl = clusts
    M = max(clusts)
    counter = 1
    for (i in 1:M) {
        if (sum(clusts == i) > 0) {
            newCl[clusts == i] = counter
            counter = counter + 1
        }
    }
    return (newCl)
}

# merge small clusters : pourrait etre interessant aussi avec du hungarian matching ;
#idea : mean of 10 smallest distances from one class to another
#complexity around K^3 (like hungarian matching...)
fusion_smcl = function(data, clusts, minszcl) {
    dists = as.matrix(dist(data))
    n = length(clusts) ; seqVect = 1:n
    newCl = clusts
    M = max(clusts) ; counter = M

    repeat {
        test = FALSE
        for (i in 1:M) {
            sI = seqVect[newCl == i]
            lI = length(sI)
            if (lI > 0 && lI < minszcl) {
                test = TRUE
                counter = counter - 1
                #too small cluster, will join another one :
                minDist = Inf ; minInd = 1
                for (j in 1:M) {
                    sJ = seqVect[newCl == j]
                    lJ = length(sJ)
                    if (j!=i && lJ > 0) {
                        srt = sort(dists[ sI,sJ ])
                        szComp = min(10, min(lI, lJ))
                        tmp = mean(srt[1:szComp])
                        if (tmp < minDist) {
                            minDist = tmp
                            minInd = j
                        }
                    }
                }
                newCl[sI] = minInd
            }
        }
        if (!test) break
    }

    return (newCl)
}

#version saying "I want exactly K clusters" from clusts :
mergeToK = function(data, clusts, K) {
    dists = as.matrix(dist(data))
    n = length(clusts) ; seqVect = 1:n
    clusts = reordering(clusts) ; M = max(clusts)

    while (M > K) { #M is decreasing one by one

        #first find smallest cluster :
        minInd = 0 ; minSize = Inf
        for (i in 1:M) {
            lI = length(seqVect[clusts == i])
            if (lI < minSize) {
                minSize = lI
                minInd = i
            }
        }
        sm = seqVect[clusts == minInd]
        lm = length(sm)
		oldMinInd=minInd

        #now merge it :
        minInd = 1 ; minDist = Inf
        for (i in 1:M) {
            if (i == oldMinInd) next
            sI = seqVect[clusts == i]
            lI = length(sI)
            srt = sort(dists[ sI,sm ])
            szComp = min(10, min(lI, lm))
            tmp = mean(srt[1:szComp])
            if (tmp < minDist) {
                minDist = tmp
                minInd = i
            }
        }

        clusts[sm] = minInd
        clusts = reordering(clusts)
        M = M-1
    }

    return (clusts)
}

