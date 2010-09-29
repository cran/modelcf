#ifndef _DIMENSION_
#define _DIMENSION_

#include <R.h>

/* listes chainees = simplexes */

typedef struct Simplex_ {
    int element;
    struct Simplex_* next;
} Simplex;

typedef struct Simplices_ {
    Simplex* current;
    struct Simplices_* next;
} Simplices;


/***************************************************/
/** "INTRUSIVE" METHODS : **/

Simplex* buildSimplex(int newElt, Simplex* simp) {
    Simplex* newSimp = (Simplex*)malloc(sizeof(Simplex));
    newSimp->element = newElt;
    newSimp->next = simp;
    return newSimp;
}

Simplex* copySimplex(Simplex* simp) {
    Simplex *tmpS, *newSimp;
    newSimp = (Simplex*)malloc(sizeof(Simplex));
    tmpS = newSimp;
    while (simp != NULL) {
        tmpS->element = simp->element;
        if (simp->next != NULL) {
            tmpS->next = (Simplex*)malloc(sizeof(Simplex));
            tmpS = tmpS->next;
        }
        else tmpS->next = NULL;
        simp = simp->next;
    }
    return newSimp;
}

Simplices* addSimplex(Simplex* simp,Simplices* simps, int duplicate) {
    Simplices* simplices;
    simplices = (Simplices*)malloc(sizeof(Simplices));
    if (duplicate) simplices->current = copySimplex(simp);
    else simplices->current = simp;
    if (simps->current == NULL) {
        simplices->next = NULL;
        free(simps); /* initialisation case */
    }
    else simplices->next = simps;
    return simplices;
}


/***************************************************/
/** "FREE MEMORY" METHODS : **/

/* free a simplex */
void freeSimplex(Simplex* simplex) {
    Simplex *tmpS;
    /* delete every elements of simplex */
    while (simplex != NULL) {
        tmpS = simplex->next;
        free(simplex);
        simplex = tmpS;
    }
}

/* free all simplices = "DESTRUCTOR" */
void freeSimplices(Simplices* simplices) {
    Simplices* tmpSS;
    /* run over simplices */
    while (simplices != NULL) {
        tmpSS = simplices->next;
        /* delete every elements of current simplex */
        freeSimplex(simplices->current);
        /* finally delete current simplex, and go to next */
        free(simplices);
        simplices = tmpSS;
    }
}


/***************************************************/
/** ESTIMATE DIMENSIONALITY : **/

/* estimate the dimensionality of a data set = MAIN FUNCTION */
void simpDim(int* neighbs, int* n, int* kNN, int* dim);

#endif
