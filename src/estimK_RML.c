#include <R.h>
#include "sorting.h"

/* estimate the number of neighbors for each data point */
void estimK_RML(double* data, int* n, int* m, int* knnmin, int* knnmax, double* tsoft, int* tabK) {

    /* declarations... */
    double *dists, *sortDists, maxDist, tmp, tmp2, dp;
    int *positions, *testNodes, *neighbNodes, i, j, l, t, countNeighb, testN;

    /* convention pour les tableaux a 2 dims : a[i*n+j] == a[i,j] */
    dists = (double*)R_alloc((*n)*(*n),sizeof(double));
    sortDists = (double*)R_alloc((*n)*(*n),sizeof(double));
    positions = (int*) R_alloc((*n)*(*n),sizeof(int));
    maxDist = 0.0;
    for (i=0; i < *n - 1; i++) {
        for (j=i+1; j < *n; j++) {
            tmp = 0.0;
            for (l=0; l < *m; l++)
                tmp += (data[i*(*m)+l] - data[j*(*m)+l]) * (data[i*(*m)+l] - data[j*(*m)+l]);
            dists[i*(*n)+j] = sqrt(tmp);
            dists[j*(*n)+i] = dists[i*(*n)+j];
            if (dists[i*(*n)+j] > maxDist) maxDist = dists[i*(*n)+j];
            sortDists[i*(*n)+j] = dists[i*(*n)+j];
            sortDists[j*(*n)+i] = dists[i*(*n)+j];
        }
    }

    /* tabK est initialise a sa taille maxi (knnmax * n) */
    neighbNodes = (int*) R_alloc(*knnmax,sizeof(int));
    testNodes = (int*) R_alloc(*knnmax,sizeof(int));

    /* initial quicksort(s) */
    for (i=0; i < *n; i++) {
        for (j=0; j < *n; j++) positions[i*(*n)+j] = j;
        sortDists[i*(*n)+i] = maxDist + 1.0;
        quickSort(sortDists,positions,i*(*n),(i+1)*(*n) - 1);
    }

    /* loops for finding neighborhoods */
    for (i=0; i < *n; i++) {
        /* reinit all nodes nearby */
        for (j=0; j < *knnmax; j++) {
            neighbNodes[j] = -1;
            testNodes[j] = -1;
        }

        countNeighb=0;
        testN = 0;
        for (j=0; j < *knnmax; j++) {
            testN = 1;
            for (l=0; l < *knnmax; l++) {
                if (testNodes[l]>=0) {
                    /* compute dot product */
                    dp = 0.0;
                    for (t=0; t < *m; t++)
                        dp += ( (data[i*(*m)+t] - data[testNodes[l]*(*m)+t]) * (data[positions[i*(*n)+j]*(*m)+t] - data[testNodes[l]*(*m)+t]) );
                    dp /= (dists[i*(*n)+testNodes[l]] * dists[positions[i*(*n)+j]*(*n)+testNodes[l]]);
                    if (dp < - *tsoft) { /* adoucissant... */
                        testN = 0;
                        break;
                    }
                }
                else break;
            }
            tabK[i*(*knnmax)+j] = 0;
            if (testN) {
                neighbNodes[countNeighb] = positions[i*(*n)+j];
                tabK[i*(*knnmax)+countNeighb] = positions[i*(*n)+j] + 1; /* j est le counterNeighb-ieme voisin de data[i,] */
                countNeighb++;
            }
            testNodes[j] = positions[i*(*n)+j];
        }

        if (countNeighb < *knnmin) {
            /* patch : il faut au moins knnmin voisins (idee : dimension..) */
            /* petit passage un peu lent : il faut trouver les plus proche voisins absents de neighbNodes */
            for (j=0; j < *knnmax; j++) {
                testN = 1;
                for (l=0; l < countNeighb; l++) {
                    if (neighbNodes[l] == positions[i*(*n)+j]) {
                        testN = 0;
                        break;
                    }
                }
                if (testN) {
                    tabK[i*(*knnmax)+countNeighb] = positions[i*(*n)+j] + 1;
                    countNeighb++;
                }
                if (countNeighb >= *knnmin) break;
            }
        }

        if (countNeighb > *knnmin) {
            /* patch : on elimine les "short-circuit" edges */
            /* a) calcul de la somme des distances */
            tmp = 0.0;
            for (j=0; j < countNeighb; j++) tmp = tmp + dists[ i*(*n) + tabK[i*(*knnmax)+j] ];

            /* b) elimination des sommets pour lesquels le ratio est >= 0.55 [0.55 assez arbitraire...] */
            for (j = countNeighb-1; j >= *knnmin; j--) {
                tmp2 = dists[ i*(*n) + tabK[i*(*knnmax)+j] ] / tmp;
                if (tmp2 >= 0.55) tabK[i*(*knnmax)+j] = 0;
                else break;
            }
        }
    }
}
