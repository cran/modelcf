#include "dimension.h"

/***************************************************/
/** "UNINTRUSIVE" METHODS : **/

/* check if a n-uple "set" has no doublon in "simplices" */
int isNotDbl(Simplex* set,Simplices* simplices) {
    int check;
    Simplex *tmpS, *tmpSet;
    if (simplices->current == NULL) return 1;
    while (simplices != NULL) {
        tmpS = simplices->current;
        tmpSet = set;
        check = 0;
        while (tmpS != NULL) {
            if (tmpS->element != tmpSet->element) {
                check = 1;
                break;
            }
            tmpS = tmpS->next;
            tmpSet = tmpSet->next;
        }
        if (check == 0) return 0;
        simplices = simplices->next;
    }
    return 1;
}

/* check if all the elements in simplex simp are in the array (of neighbors) tab */
int elemsInTab(Simplex* simp, int* tab, int firstIndTab, int lastIndTab) {
    int i,check;
    while (simp != NULL) {
        check = 0;
        for (i=firstIndTab; i<lastIndTab; i++) {
            if (simp->element == tab[i]) {
                check = 1;
                break;
            }
        }
        if (!check) return 0;
        simp = simp->next;
    }
    return 1;
}


/***************************************************/
/** ESTIMATE DIMENSIONALITY : **/

/* estimate the dimensionality of a data set = MAIN FUNCTION */
void simpDim(int* neighbs, int* n, int* kNN, int* dim) {
    int i, j, curDim, countS;
    Simplices *simplices, *newSimplices, *tmpSS;
    Simplex *tmpS1, *tmpS2;

    /* first step : find all the 1D simplices (2 elements) */
    simplices = (Simplices*)malloc(sizeof(Simplices));
    simplices->current = NULL;
    simplices->next = NULL;
    for (i=0; i < *n; i++) {
        for (j=0; j < *kNN; j++) {
            if (neighbs[i * (*kNN) + j] > 0) {
                tmpS1 = (Simplex*)malloc(sizeof(Simplex));
                tmpS1->element = neighbs[i * (*kNN) + j];
                tmpS2 = (Simplex*)malloc(sizeof(Simplex));
                tmpS2->element = i+1;
                tmpS1->next = tmpS2;
                tmpS2->next = NULL;

                /* on ajoute ce simplexe (2 elements) */
                if (isNotDbl(tmpS1,simplices))
                    simplices = addSimplex(tmpS1,simplices,0);
                else {
                    free(tmpS1);
                    free(tmpS2);
                }
            }
            else break;
        }
    }

    /* second step : find all the higher-order simplices iteratively */
    curDim = 2;
    while (1) {
        newSimplices = (Simplices*)malloc(sizeof(Simplices));
        newSimplices->current = NULL;
        newSimplices->next = NULL;

        /* search for (curDim+1)-dimensional simplices ; print stats first */
        countS=0;
        tmpSS = simplices;
        while (tmpSS != NULL) {
            countS = countS + 1;
            tmpSS = tmpSS->next;
        }
        tmpSS = simplices;
        printf("%i simplices at dimension %i\n",countS,curDim-1);

        while (tmpSS != NULL) {
            /* on parcourt tous les elements 1--> n en cherchant si ceux de tmpSS->current sont dans les voisins */
            for (i=0; i < *n; i++) {

                /* cas 1 : on y arrive pas ; on ne fait donc rien */
                if (elemsInTab(tmpSS->current, neighbs, i*(*kNN), (i+1)*(*kNN))) {
                    tmpS1 = buildSimplex(i+1,tmpSS->current);

                    /* cas 2 : on y arrive, mais c'est un doublon dans newSimplices --> on ne fait rien */
                    if (isNotDbl(tmpS1,newSimplices)) {

                        /* cas 3: on y arrive et c'est pas un doublon : on l'ajoute a newSimplices */
                        newSimplices = addSimplex(tmpS1,newSimplices,1);
                    }
                    free(tmpS1);
                }
            }
            tmpSS = tmpSS->next;
        }

        /* free the old simplices */
        freeSimplices(simplices);

        /* if no simplices have been found, the estimated dimension is curDim - 1 */
        if (newSimplices->current == NULL) {
            *dim = curDim-1;
            free(newSimplices);
            return;
        }
        curDim = curDim + 1;
        simplices = newSimplices;
    }
    return;
}
