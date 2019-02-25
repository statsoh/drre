#include <R.h>
#include "nrutilR.h"
#include <math.h>

#define FREERETURN2 {free_dvector(objtemptemp,0,P);return objtemp;}
double objfun(int i, double *beta, int n, int M, int P,int MCol, double *u, double *dt, int *MNum,double **z, int **MMat, double **B){

    int r,g,gg;
    double lin=0.0,lin2=0.0,objtemp=0.0,elin=0.0, *objtemptemp;
    objtemptemp=dvector(0,P);
    for (r=0;r<=P;r++) objtemptemp[r]=0;


///Get: lin, exp(lin), then objtemp

    if (i==0){ /// objective function for i==0
        for (r=0;r<n;r++){
            if (u[i] < dt[r]){  /// dt = cen<---- censoring time, not last event time.
                lin=0.0;
                for (g=0;g<=P;g++) lin+= z[r][g]*beta[g];
                elin=exp(lin);
                objtemp-=elin;
            }
        }
        for (g=0;g<=P;g++){
            for (gg=0;gg<MNum[i];gg++){
                objtemp+=z[MMat[i][gg]-1][g]*beta[g];
            }
        }
        
    }else{ /// objective function for i!=0
        for (r=0;r<n;r++){
            if (u[i] < dt[r]){
                lin=0;
                lin2=0;
                for (g=0;g<=P;g++){
                    lin+= z[r][g]*B[i][g];
                    lin2+= z[r][g]*beta[g];
                }
                elin=exp(lin);
                objtemp-=exp(lin2);
                for (g=0;g<=P;g++){
                    objtemptemp[g]+=z[r][g]*elin;
                }
            }
        }
        for (g=0;g<=P;g++){
            objtemp+= objtemptemp[g]*beta[g];
        }
        for (g=0;g<=P;g++){
            for (gg=0;gg<MNum[i];gg++){
                objtemp+=z[MMat[i][gg]-1][g]*beta[g];
            }
        }
    }
        /// Objective function is convex. Maximum is where the root exists.

    FREERETURN2;
}
