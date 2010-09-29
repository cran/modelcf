/********************************************************************
 ********************************************************************
 **
 ** libhungarian by Cyrill Stachniss, 2004
 ** -- modified (very) slightly by Benjamin Auder, 2010
 ** -- (verbose and printing routines deleted, because used in a R package)
 **
 **
 ** Solving the Minimum Assignment Problem using the 
 ** Hungarian Method.
 **
 ** ** This file may be freely copied and distributed! **
 **
 ** Parts of the used code was originally provided by the 
 ** "Stanford GraphGase", but I made changes to this code.
 ** As asked by  the copyright node of the "Stanford GraphGase", 
 ** I hereby proclaim that this file are *NOT* part of the
 ** "Stanford GraphGase" distribution!
 **
 ** This file is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied 
 ** warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 ** PURPOSE.
 **
 ********************************************************************
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "hungarian.h"

#define INF (0x7FFFFFFF)
#define hungarian_test_alloc(X) do {if ((void *)(X) == NULL) fprintf(stderr, "Out of memory in %s, (%s, line %d).\n", __FUNCTION__, __FILE__, __LINE__); } while (0)

void hungarian_init(hungarian_problem_t* p, int* cost_matrix, int rows, int cols, int mode) {

    int i,j, org_cols, org_rows;
    int max_cost;
    max_cost = 0;
    
    org_cols = cols;
    org_rows = rows;
    
    // is the number of cols  not equal to number of rows ? 
    // if yes, expand with 0-cols / 0-cols
    rows = (cols < rows) ? rows : cols;
    cols = rows;
    
    p->num_rows = rows;
    p->num_cols = cols;
    
    p->cost = (int**)calloc(rows,sizeof(int*));
    hungarian_test_alloc(p->cost);
    p->assignment = (int**)calloc(rows,sizeof(int*));
    hungarian_test_alloc(p->assignment);

    for(i=0; i<p->num_rows; i++) {
        p->cost[i] = (int*)calloc(cols,sizeof(int));
        hungarian_test_alloc(p->cost[i]);
        p->assignment[i] = (int*)calloc(cols,sizeof(int));
        hungarian_test_alloc(p->assignment[i]);
        for(j=0; j<p->num_cols; j++) {
        p->cost[i][j] =  (i < org_rows && j < org_cols) ? cost_matrix[i*(p->num_cols)+j] : 0;
        p->assignment[i][j] = 0;
    
        if (max_cost < p->cost[i][j])
        max_cost = p->cost[i][j];
        }
    }

    if (mode == HUNGARIAN_MODE_MAXIMIZE_UTIL) {
        for(i=0; i<p->num_rows; i++) {
            for(j=0; j<p->num_cols; j++) {
                p->cost[i][j] =  max_cost - p->cost[i][j];
            }
        }
    }
}

void hungarian_free(hungarian_problem_t* p) {
    int i;
    for(i=0; i<p->num_rows; i++) {
        free(p->cost[i]);
        free(p->assignment[i]);
    }
    free(p->cost);
    free(p->assignment);
    p->cost = NULL;
    p->assignment = NULL;
}

void hungarian_solve(hungarian_problem_t* p) {

    int i, j, m, n, k, l, s, t, q, unmatched, cost;
    int* col_mate;
    int* row_mate;
    int* parent_row;
    int* unchosen_row;
    int* row_dec;
    int* col_inc;
    int* slack;
    int* slack_row;
    
    cost=0;
    m =p->num_rows;
    n =p->num_cols;
    
    col_mate = (int*)calloc(p->num_rows,sizeof(int));
    hungarian_test_alloc(col_mate);
    unchosen_row = (int*)calloc(p->num_rows,sizeof(int));
    hungarian_test_alloc(unchosen_row);
    row_dec  = (int*)calloc(p->num_rows,sizeof(int));
    hungarian_test_alloc(row_dec);
    slack_row  = (int*)calloc(p->num_rows,sizeof(int));
    hungarian_test_alloc(slack_row);
    
    row_mate = (int*)calloc(p->num_cols,sizeof(int));
    hungarian_test_alloc(row_mate);
    parent_row = (int*)calloc(p->num_cols,sizeof(int));
    hungarian_test_alloc(parent_row);
    col_inc = (int*)calloc(p->num_cols,sizeof(int));
    hungarian_test_alloc(col_inc);
    slack = (int*)calloc(p->num_cols,sizeof(int));
    hungarian_test_alloc(slack);

    for (i=0;i<p->num_rows;i++) {
        col_mate[i]=0;
        unchosen_row[i]=0;
        row_dec[i]=0;
        slack_row[i]=0;
    }
    for (j=0;j<p->num_cols;j++) {
        row_mate[j]=0;
        parent_row[j] = 0;
        col_inc[j]=0;
        slack[j]=0;
    }

    for (i=0;i<p->num_rows;++i)
        for (j=0;j<p->num_cols;++j)
            p->assignment[i][j]=HUNGARIAN_NOT_ASSIGNED;

    // Begin subtract column minima in order to start with lots of zeroes 12
    for (l=0;l<n;l++) {
        s=p->cost[0][l];
        for (k=1;k<m;k++) 
        if (p->cost[k][l]<s)
        s=p->cost[k][l];
        cost+=s;
        if (s!=0)
        for (k=0;k<m;k++)
        p->cost[k][l]-=s;
    }
    // End subtract column minima in order to start with lots of zeroes 12

    // Begin initial state 16
    t=0;
    for (l=0;l<n;l++) {
        row_mate[l]= -1;
        parent_row[l]= -1;
        col_inc[l]=0;
        slack[l]=INF;
    }
    for (k=0;k<m;k++) {
        s=p->cost[k][0];
        for (l=1;l<n;l++)
        if (p->cost[k][l]<s)
        s=p->cost[k][l];
        row_dec[k]=s;
        for (l=0;l<n;l++)
        if (s==p->cost[k][l] && row_mate[l]<0) {
            col_mate[k]=l;
            row_mate[l]=k;
            goto row_done;
        }
        col_mate[k]= -1;
        unchosen_row[t++]=k;
        row_done:
        ;
    }
    // End initial state 16

    // Begin Hungarian algorithm 18
    if (t==0)
        goto done;
    unmatched=t;
    while (1) {
        q=0;
        while (1) {
            while (q<t) {
                // Begin explore node q of the forest 19
                {
                k=unchosen_row[q];
                s=row_dec[k];
                for (l=0;l<n;l++)
                    if (slack[l]) {
                        int del;
                        del=p->cost[k][l]-s+col_inc[l];
                        if (del<slack[l]) {
                            if (del==0) {
                                if (row_mate[l]<0)
                                    goto breakthru;
                                slack[l]=0;
                                parent_row[l]=k;
                                unchosen_row[t++]=row_mate[l];
                            }
                            else {
                                slack[l]=del;
                                slack_row[l]=k;
                            }
                        }
                    }
                }
                // End explore node q of the forest 19
                q++;
            }

            // Begin introduce a new zero into the matrix 21
            s=INF;
            for (l=0;l<n;l++)
                if (slack[l] && slack[l]<s)
                    s=slack[l];
            for (q=0;q<t;q++)
                row_dec[unchosen_row[q]]+=s;
            for (l=0;l<n;l++) {
                if (slack[l]) {
                    slack[l]-=s;
                    if (slack[l]==0) {
                        // Begin look at a new zero 22
                        k=slack_row[l];
                        if (row_mate[l]<0) {
                            for (j=l+1;j<n;j++)
                                if (slack[j]==0)
                                    col_inc[j]+=s;
                            goto breakthru;
                        }
                        else {
                            parent_row[l]=k;
                            unchosen_row[t++]=row_mate[l];
                        }
                        // End look at a new zero 22
                    }
                }
                else
                    col_inc[l]+=s;
                    // End introduce a new zero into the matrix 21
            }
        }
        breakthru:
        // Begin update the matching 20
        while (1) {
            j=col_mate[k];
            col_mate[k]=l;
            row_mate[l]=k;
            if (j<0)
            break;
            k=parent_row[j];
            l=j;
        }
        // End update the matching 20
        if (--unmatched==0)
            goto done;
        // Begin get ready for another stage 17
        t=0;
        for (l=0;l<n;l++) {
            parent_row[l]= -1;
            slack[l]=INF;
        }
        for (k=0;k<m;k++)
            if (col_mate[k]<0)
                unchosen_row[t++]=k;
        // End get ready for another stage 17
    }

    done:

    // Begin doublecheck the solution 23
    for (k=0;k<m;k++)
        for (l=0;l<n;l++)
        if (p->cost[k][l]<row_dec[k]-col_inc[l])
        exit(0);
    for (k=0;k<m;k++) {
        l=col_mate[k];
        if (l<0 || p->cost[k][l]!=row_dec[k]-col_inc[l])
        exit(0);
    }
    k=0;
    for (l=0;l<n;l++)
        if (col_inc[l])
            k++;
    if (k>m)
        exit(0);
    // End doublecheck the solution 23
    // End Hungarian algorithm 18

    for (i=0;i<m;++i) {
        p->assignment[i][col_mate[i]]=HUNGARIAN_ASSIGNED;
        /*TRACE("%d - %d\n", i, col_mate[i]);*/
    }
    for (k=0;k<m;++k) {
        for (l=0;l<n;++l) {
            p->cost[k][l]=p->cost[k][l]-row_dec[k]+col_inc[l];
        }
    }
    for (i=0;i<m;i++)
        cost+=row_dec[i];
    for (i=0;i<n;i++)
        cost-=col_inc[i];

    free(slack);
    free(col_inc);
    free(parent_row);
    free(row_mate);
    free(slack_row);
    free(row_dec);
    free(unchosen_row);
    free(col_mate);
}

/* get the optimal assignment, by calling hungarian_solve above */
void getOptAssign(int* P1, int* P2, int* maxInd, int* n, int* assign_cl) {

    /* first, determine weights by computing intersections
     * reads like: weights[i*(*n)+j] == gain when cluster i in P1
     * is assigned to cluster j in P2
     */
    int i,j,k,inter;
    int* utils = malloc((*maxInd) * (*maxInd) * sizeof(int));
    for (i=1; i <= *maxInd; i++) {
        for (j=1; j <= *maxInd; j++) {
            //try label i with label j
            inter=0;
            for (k=0; k < *n; k++) {
                if (P1[k]==i && P2[k]==j) inter++;
            }
            utils[(i-1)*(*maxInd)+(j-1)] = inter;
        }
    }

    //then solve problem :
    hungarian_problem_t p;
    hungarian_init(&p, utils, *maxInd, *maxInd, HUNGARIAN_MODE_MAXIMIZE_UTIL) ;
    hungarian_solve(&p);

    //and now get optimal assignment :
    for (i=0; i < *maxInd; i++) {
        for (j=0; j < *maxInd; j++) {
            if (p.assignment[i][j]) assign_cl[i] = j;
        }
    }
    hungarian_free(&p);
}

