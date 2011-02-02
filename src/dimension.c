#include <R.h>
#include <stdlib.h>
#include <time.h>

//#include <stdio.h>

/* "random" hierarchical clustering : */
void rand_hclust(double* dists, int* n, double* infini, int* stats) {
    int i, j, I, J, k, l, s, t, n2 = (*n)*(*n);
    int* nonempty = (int*) malloc((*n) * sizeof(int));
    for (i=0; i < *n; i++) nonempty[i] = 1;

    int* classes = (int*)calloc(n2, sizeof(int));
    /* initially, each element is a cluster : */
    for (i=0; i < *n; i++) classes[i*(*n)] = i;
    int* indxNext = (int*) malloc((*n) * sizeof(int));
    for (i=0; i < *n; i++) indxNext[i] = 1;

    double* dissims = (double*) malloc(n2 * sizeof(double));
    for (i=0; i<n2; i++) dissims[i] = dists[i];
    double* candsFusion = (double*) malloc(n2 * sizeof(double));
    int isInfty, nsize, tmp, nbCouples, single, chosen;
    int* couples = (int*) malloc( (*n) * (*n - 1) * sizeof(int) );
    int* couple = (int*) malloc(2 * sizeof(int));
    srand ( time(NULL) );
    double M;

    for (nsize = *n; nsize>0; nsize--) {
        for (i=0; i<n2; i++) candsFusion[i] = *infini;
        int* nonvals = (int*) malloc(nsize * sizeof(int)); // = refVect[nonempty]
        tmp = 0;
        for (i=0; i < *n; i++) {
            if (nonempty[i]) {
                nonvals[tmp] = i;
                tmp++;
            }
        }

        /* recherche des candidats à la fusion */
        M = 0.0;
        for (i=0; i<nsize; i++) {
            for (j=0; j<nsize; j++) {
                // candsFusion[I,J] = 1.0 / mean( dissims[ classes[[j]],classes[[k]] ] )
                I = nonvals[i];
                J = nonvals[j];
                isInfty = 0;
                if (I < J) {
                    candsFusion[I*(*n)+J] = 0.0;
                    for (k=0; k<indxNext[I]; k++) {
                        for (l=0; l<indxNext[J]; l++) {
                            if (dissims[ classes[I*(*n)+k]*(*n) + classes[J*(*n)+l] ] == *infini) {
                                isInfty = 1;
                                break;
                            }
                            candsFusion[I*(*n)+J] += dissims[ classes[I*(*n)+k]*(*n) + classes[J*(*n)+l] ];
                        }
                        if (isInfty) break;
                    }
                    if (isInfty) candsFusion[I*(*n)+J] = 0.0;
                    else {
                        candsFusion[I*(*n)+J] /= (indxNext[I] * indxNext[J]);
                        candsFusion[I*(*n)+J] = 1.0/candsFusion[I*(*n)+J];
                        if (candsFusion[I*(*n)+J] > M) M = candsFusion[I*(*n)+J];
                    }
                }
            }
        }

        if (M == 0.0) {
            for (i=0; i < *n; i++) {
                if (indxNext[i] >= 1) stats[indxNext[i]-1]++;
            }
            free(nonempty);
            free(classes);
            free(indxNext);
            free(dissims);
            free(candsFusion);
            free(couples);
            free(couple);
            free(nonvals);
            return;
        }

        /* algo rejet pour picker un couple parmi les éligibles : */
        nbCouples = 0;
        for (i=0; i<nsize; i++) {
            for (j=0; j<nsize; j++) {
                I = nonvals[i];
                J = nonvals[j];
                if (I < J && candsFusion[I*(*n)+J] > 0.0) {
                    couples[2*nbCouples] = I;
                    couples[2*nbCouples+1] = J;
                    nbCouples++;
                }
            }
        }
        free(nonvals);

        // on constitue une liste de couples, puis on en tire un au hasard (uniforme)
        // puis on tire VA U[0,M], et si <= "proba" du couple (celle dans candsFusion) alors OK
        single = 1;
        chosen = -1;
        while (single) {
            chosen = rand() % nbCouples;
            if (M*(double)rand()/RAND_MAX <= candsFusion[ couples[2*chosen]*(*n) + couples[2*chosen+1] ])
                single = 0;
        }
        couple[0] = couples[2*chosen];
        couple[1] = couples[2*chosen+1];

        /* Puis fusion (mettre à jour nonempty, classes et dissims) */
        s = couple[0];
        t = couple[1];
        nonempty[t] = 0;
        for (i=0; i<indxNext[t]; i++) classes[s*(*n)+indxNext[s]+i] = classes[t*(*n)+i];
        indxNext[s] += indxNext[t];
        indxNext[t] = 0;
        for (i=0; i < *n; i++) {
            if (dissims[s*(*n)+i] == *infini || dissims[t*(*n)+i] == *infini)
                dissims[s*(*n)+i] = *infini;
            else {
                if (dissims[s*(*n)+i] > dissims[t*(*n)+i])
                    dissims[s*(*n)+i] = dissims[t*(*n)+i];
            }
            dissims[i*(*n)+s] = dissims[s*(*n)+i];
        }
    }

    free(nonempty);
    free(classes);
    free(indxNext);
    free(dissims);
    free(candsFusion);
    free(couples);
    free(couple);

    return;
}

