#include <R.h>

/* check if a vector contain any 0 */
int checkZero(int* vect, int n) {
    int i;
    for (i=0; i<n; i++) {
        if (vect[i] == 0) return i;
    }
    return 0;
}

/* connex components analysis for WEAK definition, mutual-kNN graph */
void weakConnexSym(double* binMat, int n, int* cc) {
    int* alreadyExpanded;
    int curInd, label, hasChanged, hasZeros, i, j;

    for (i=0; i < n; i++) cc[i] = 0;
    curInd = 0;
    hasZeros = 1;
    alreadyExpanded = (int*) R_alloc (n, sizeof(int));

    /* while the entire graph hasn't been visited */
    while (hasZeros) {
        hasChanged = 1;
        label = curInd+1;
        cc[curInd] = label;
        for (i=curInd; i < n; i++) alreadyExpanded[i] = 0;

        /* while new elements are discovered in current component */
        while (hasChanged) {
            hasChanged = 0;
            for (i=curInd; i < n; i++) {
                if (cc[i] == label && !alreadyExpanded[i]) {
                    for (j=curInd+1; j < n; j++) {
                        if (!cc[j] && binMat[i*n+j]) {
                            hasChanged=1;
                            cc[j] = label;
                        }
                    }
                    alreadyExpanded[i] = 1;
                }
            }
        }

        /* curInd is set to the next unexplored index (if possible) */
        curInd = checkZero(cc,n);
        hasZeros = curInd;
    }
}

/* aux function for "assymmetric" connexity : expands local neighborhoods */
void expand(double* binMat, int n, int* cc, int index) {
    int i, minLab;

    /* find all the direct neighbors of i and compute their min. label */
    minLab = cc[index];
    for (i=0; i < n; i++) {
        if (binMat[index*n + i] && cc[i] < minLab)
            minLab = cc[i];
    }

    /* assign the lowest label to all of them */
    for (i=0; i < n; i++) {
        if (binMat[index*n + i]) cc[i] = minLab;
    }
}

/* connex components analysis for WEAK definition, "assymetric" graph */
void weakConnexAssym(double* binMat, int n, int* cc) {
    int* oldCC = (int*) R_alloc (n, sizeof(int));
    int i, ccHasChanged;

    /* strategy = local neighbors expanding */
    for (i=0; i < n; i++) cc[i] = i+1;

    ccHasChanged = 1;
    while (ccHasChanged) {

        oldCC = cc;
        for (i=0; i < n; i++) expand(binMat, n, cc, i);
        ccHasChanged = 0;
        for (i=0; i < n; i++) {
            if (cc[i] != oldCC[i]) {
                ccHasChanged = 1;
                break;
            }
        }
    }
}

/* explore the connectivity of a graph (binarySim = adjacency matrix) */
void getConnex(double* binMat, int* n, int* ctype, int* cc) {
    if (*ctype) weakConnexSym(binMat, *n, cc);
    else weakConnexAssym(binMat, *n, cc);
}

