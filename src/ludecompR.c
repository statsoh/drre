#include <R.h>
#include <math.h>
#include "nrutilR.h"
#define TINY 1.0e-20;   // A small number


void ludcmp(double **a, int dim, int *indx, double *d){
// Given a matrix a[1..dim][1..dim], this routine replaces it by the LU decomposition of a rowwise permutation of itself.
// a and dim are input.
// a is output, arranged as in equation (2.3.14);
// indx[1..dim] is an output vector that records the row permutation effected by the partial pivoting
// d is output as +-1 depending on whether the number of row interchanges was even or odd, respectively.
// This routine is used in compination with lubksb to solve linear equations or invert a matrix.

    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;      // vv stores the implicit scaling of each row.

    vv=dvector(1,dim);
    *d=1.0;         // no row interchanges yet
    for (i=1;i<=dim;i++){
                    // Loop over rows to get the implicit scaling information.
        big=0.0;
        for (j=1;j<=dim;j++)
            if ((temp=fabs(a[i][j]))> big) big=temp;
        if (big==0.0) warning("Singular matrix in routine LUDCMP"); // no nonzero largest element.
        vv[i]=1.0/big;      // save the scaling
    }
    /////////////// This is the loop over columns of Crout's method.
    for (j=1;j<=dim;j++){
        for (i=1;i<j;i++){  // This is equation (2.3.12) except for i=j!!!!!!!!!!!!!!
            sum=a[i][j];
            for (k=1;k<i;k++) sum-=a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;            // Initialize for the search for largest pivot element
        for (i=j;i<=dim;i++){ // This is i=j of equation (2.3.12) and i=j+1...N of eq (2.3.13)
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum-=a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big){
                // Is the figure of merit for the pivot better than the best so far?
                big=dum;
                imax=i;
            }
        }
        if (j!= imax){   // Do we need to interchange rows?
            for(k=1;k<=dim;k++){  // Yes, do so...
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d= -(*d);          // and change the parity of d
            vv[imax]=vv[j];     // Also interchange the scale factor.
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;   // If the pivot element is zero the matrix is singular
        if (j!= dim){
            // Now finally, divide by the pivot element.
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=dim;i++) a[i][j] *= dum;
        }
    }
    free_dvector(vv,1,dim);
}

void lubksb(double **a, int dim, int *indx, double *b){
// Solve the set of 'dim' linear equations AX=B. Here a[1..dim][1..dim] is input, not as the matrix A but rather as its LU decomposition, determined by the routine ludcmp.
// indx[1..dim] is input as the permutation vector returned by ludcomp.
// b[1..dim] is input as the right-hand side vector B, and returns with the solution vector X.
// a, dim, and indx are NOT modified by this routine and can be left in place for successive calls with different right-hand side b.

    int i,ii=0,ip,j;
    double sum;

    for (i=1;i<=dim;i++){
        // When ii is set to a positive value, it will become the index of the first nonvanishing element of b.
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
    }
    for (i=dim;i>=1;i--){
        sum=b[i];
        for (j=i+1;j<=dim;j++) sum-= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}
