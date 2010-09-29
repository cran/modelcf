#include "sorting.h"

/* fonction auxiliaire pour le quicksort */
int partitionner(double *tableau, int* positions, int p, int r) {
    double pivot = tableau[p], tmpF;
    int i = p-1, j = r+1, tmpI;

    while (1) {
        do j--;
        while (tableau[j] > pivot);

        do i++;
        while (tableau[i] < pivot);

        if (i < j) {
            tmpF = tableau[i];
            tableau[i] = tableau[j];
            tableau[j] = tmpF;
            tmpI = positions[i];
            positions[i] = positions[j];
            positions[j] = tmpI;
        }
        else return j;
    }
    return j;
}

void quickSort(double *tableau, int* positions, int p, int r) {
    int q;
    if (p < r) {
        q = partitionner(tableau, positions, p, r);
        quickSort(tableau, positions, p, q);
        quickSort(tableau, positions, q+1, r);
    }
}
