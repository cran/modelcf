#ifndef _SORTING_
#define _SORTING_

#include <R.h>

/* fonction auxiliaire pour le quicksort */
int partitionner(double *tableau, int* positions, int p, int r);

void quickSort(double *tableau, int* positions, int p, int r);

#endif
