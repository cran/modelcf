#auxiliar function for the equation to solve (see Lin and Zha article)
psi_RML = function(x, a2, b2, s2) {
    s2Px = s2 + x
    vB = b2 / s2Px
    vS = s2 / s2Px
    return (sum(vB * vS) - a2)
}

#auxiliar function for parameters optimization
evalParams = function(x, a2, Am, bv) {
    y = ( sqrt(a2) / sqrt(sum(x^2)) ) * x
    return (sum( (Am%*%y - bv)^2 ))
}

#solver : find the coordinates of data[curInd,] given the known embeddings
#predInd = index of the predecessor of curInd on a shortest path.
solveConstraints = function(data, embedding, curInd, predInd, predNeighbs, dim) {

    normQR = sqrt(sum((data[curInd,] - data[predInd,])^2))
    x = 0
    #get A and b (see article)
    if (dim > 1) {
        A = t( t(embedding[predNeighbs,]) - embedding[predInd,] )		
        A = A / (normQR * sqrt(rowSums ( ( t( t(embedding[predNeighbs,]) - embedding[predInd,] ) )^2) ))
        b = colSums( (data[curInd,] - data[predInd,]) * ( t(data[predNeighbs,]) - data[predInd,] ) )
        b = b / (normQR * (sqrt( rowSums( t( ( t(data[predNeighbs,]) - data[predInd,])^2 ) ) )) )
        alpha2 = normQR^2

        #if A is well-conditioned, solve equation with uniroot :
        s=svd(A,nv=0)
        if (min(s$d) > EPS() ) {
            sigma2 = (s$d)^2
            beta2 = ( t(s$u) %*% b )^2
            Ms2 = max(sigma2)
            upBound = Ms2
            while (psi_RML(upBound, alpha2, beta2, sigma2) > EPS() ) upBound = upBound + Ms2
            solveEq = uniroot( psi_RML, c(-Ms2 + EPS(), upBound), a2=alpha2, b2=beta2, s2=sigma2 )
            lambda = solveEq$root

            #now we can get x :
            I = diag(1,dim)
            tmpMat = mppsinv( t(A) %*% A + lambda * I )
            x = tmpMat %*% t(A) %*% b
        }
        else {
            #direct numerical optimization :
            op = optim(rep(1.0,dim), evalParams, a2=alpha2, Am=A, bv=b)
            x = ( sqrt(alpha2) / sqrt(sum( (op$par)^2 )) ) * op$par
        }
    }
    else {
        #special case dim==1 :
        A = embedding[predNeighbs,] - embedding[predInd,]
        A = A / (normQR * abs(embedding[predNeighbs,] - embedding[predInd,]) )
        b = colSums( (data[curInd,] - data[predInd,]) * ( t(data[predNeighbs,]) - data[predInd,] ) )
        b = b / (normQR * (sqrt( rowSums( t( ( t(data[predNeighbs,]) - data[predInd,])^2 ) ) )) )

        #shortcut to parameter estimation in 1-dim case :
        tmp1 = sum(abs(normQR*A - b)) ; tmp2 = sum(abs( - normQR*A - b))
        x=0.0
        if (tmp1 < tmp2) x = normQR
        else x = -normQR
    }
    return (x + embedding[predInd,])
}

