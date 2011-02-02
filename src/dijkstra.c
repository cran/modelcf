#include <R.h>
#include <stdlib.h>

/* Dijkstra from index debut ; dist au format vecteur ; "output" geodesics and predecessors in geodsPreds */
void dijkstra(double* dists, int* n, int* start, double* geodsPredsLevels) {

    int n1, n2, i, ind_n12;
    int* unseen;
    double infDist, minGeods;

    /* initalisations I : */
    unseen = (int*) malloc ((*n) * sizeof(int));
    infDist = 0.0;
    for (i=0; i < (*n)*(*n); i++) {
        if (dists[i] > infDist) infDist = dists[i];
    }
    infDist = (*n)*infDist + 1.0;

    /* initialisations II */
    for (i=0; i < *n; i++) {
        geodsPredsLevels[i] = infDist; /* "infini" */
        unseen[i] = 1; /* nothing seen so far */
        geodsPredsLevels[i + *n] = 0;
        geodsPredsLevels[i + 2*(*n)] = 0;
    }
    geodsPredsLevels[*start - 1] = 0.0;
    n1 = 0; /* avoid compiler warning */

    while (1) {

        /* n1 <-- le noeud dans "PasEncoreVu" avec distance geodesique la plus petite */
        minGeods = infDist;
        for (i=0; i < *n; i++) {
            if (unseen[i]) {
                if (geodsPredsLevels[i] < minGeods) {
                    n1 = i;
                    minGeods = geodsPredsLevels[i];
                }
            }
        }

        if (minGeods == infDist) break; /* all is explored */
        /* on enleve n1 de pasEncoreVus : */
        unseen[n1] = 0;

        /* Pour n2 parcourant fils(n1) (Les noeuds relies a n1 par un arc) : */
        for (n2 = 0; n2 < *n; n2++) {
            ind_n12 = n1*(*n)+n2;
            if (dists[ind_n12] > 0.0) {
                if (geodsPredsLevels[n2] > geodsPredsLevels[n1] + dists[ind_n12]) { /* distance correspond au poids de l'arc reliant n1 et n2 */
                    geodsPredsLevels[n2] = geodsPredsLevels[n1] + dists[ind_n12];
                    geodsPredsLevels[n2 + *n] = n1+1; /* Dit que pour aller a n2, il faut passer par n1 */
                    geodsPredsLevels[n2 + 2*(*n)] = geodsPredsLevels[n1 + 2*(*n)] + 1;
                }
            }
        }
    }
    free(unseen);
}

