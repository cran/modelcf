#include <R.h>
#include <math.h>
#include <stdlib.h>
#include "sorting.h"

#include <stdio.h>

/* estimate the number of neighbors for each data point */
void estimK_RML(double* data, double* rgdists, int* n, int* m, int* kmin, int* kmax, double* tsoft, int* tabK) {

    /* declarations... */
    double *dists, *sortDists, maxDist, tmp, tmp2, dp, n2;
    int *positions, *testNodes, *neighbNodes, i, j, l, t, countNeighb, testN;

    /* convention pour les tableaux a 2 dims : a[i*n|m+j] == a[i,j] */
    n2 = (*n) * (*n);
    dists = (double*)malloc(n2*sizeof(double));
    sortDists = (double*)malloc(n2*sizeof(double));
    positions = (int*) malloc(n2*sizeof(int));
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

    /* tabK est initialise a sa taille maxi (kmax * n) */
    neighbNodes = (int*) malloc((*kmax) * sizeof(int));
    testNodes = (int*) malloc((*kmax) * sizeof(int));

    /* initial quicksort(s) */
    for (i=0; i < *n; i++) {
        for (j=0; j < *n; j++) positions[i*(*n)+j] = j;
        sortDists[i*(*n)+i] = maxDist + 1.0;
        quickSort(sortDists,positions,i*(*n),(i+1)*(*n) - 1);
    }

    /* loops for finding neighborhoods */
    for (i=0; i < *n; i++) {
        /* reinit all nodes nearby */
        for (j=0; j < *kmax; j++) {
            neighbNodes[j] = -1;
            testNodes[j] = -1;
        }

        countNeighb=0;
        testN = 0;
        for (j=0; j < *kmax; j++) {
            testN = 1;
            for (l=0; l < *kmax; l++) {
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
            tabK[i*(*kmax)+j] = 0;
            if (testN) {
                neighbNodes[countNeighb] = positions[i*(*n)+j];
                tabK[i*(*kmax)+countNeighb] = positions[i*(*n)+j] + 1; /* j est le counterNeighb-ieme voisin de data[i,] */
                countNeighb++;
            }
            testNodes[j] = positions[i*(*n)+j];
        }

        if (countNeighb < *kmin) {
            /* patch : il faut au moins kmin voisins (idee : dimension..) */
            /* petit passage un peu lent : il faut trouver les plus proche voisins absents de neighbNodes */
            for (j=0; j < *kmax; j++) {
                testN = 1;
                for (l=0; l < countNeighb; l++) {
                    if (neighbNodes[l] == positions[i*(*n)+j]) {
                        testN = 0;
                        break;
                    }
                }
                if (testN) {
                    tabK[i*(*kmax)+countNeighb] = positions[i*(*n)+j] + 1;
                    countNeighb++;
                }
                if (countNeighb >= *kmin) break;
            }
        }

		if (rgdists[0] >= 0.0) {
			/* Advanced heuristic : */
			
			if (countNeighb > *kmin) {
				/* patch : on elimine les "short-circuit" edges */
				/* a) calcul de la moyenne des distances geodesiques aux plus proches voisins "de visibilite" */
				tmp = 0.0;
				for (j=0; j < countNeighb; j++) tmp = tmp + rgdists[ i*(*n) + tabK[i*(*kmax)+j] - 1 ];
				tmp /= countNeighb;

				/* b) elimination des sommets pour lesquels le ratio a la moyenne est > 1.8 [1.8 assez arbitraire...] */
				for (j = countNeighb-1; j >= *kmin; j--) {
					tmp2 = rgdists[ i*(*n) + tabK[i*(*kmax)+j] - 1 ] / tmp;
					if (tmp2 >= 1.8) tabK[i*(*kmax)+j] = 0;
					else break;
				}
			}
		}
		
		else {
			/* Basic heuristic : */
			
			if (countNeighb > *kmin) {
				/* patch : on elimine les "short-circuit" edges */
				/* a) calcul de la somme des distances geodesiques aux plus proches voisins "euclidiens" */
				tmp = 0.0;
				for (j=0; j < countNeighb; j++) tmp = tmp + dists[ i*(*n) + tabK[i*(*kmax)+j] - 1 ];
				tmp /= countNeighb;

				/* b) elimination des sommets pour lesquels le ratio a la moyenne est > 1.5 [1.5 assez arbitraire...] */
				for (j = countNeighb-1; j >= *kmin; j--) {
					tmp2 = dists[ i*(*n) + tabK[i*(*kmax)+j] - 1 ] / tmp;
					if (tmp2 >= 1.5) tabK[i*(*kmax)+j] = 0;
					else break;
				}
			}
		}
    }

    free(dists);
    free(sortDists);
    free(positions);
    free(testNodes);
    free(neighbNodes);
}

