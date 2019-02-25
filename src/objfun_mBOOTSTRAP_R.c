#include <R.h>
#include "nrutilR.h"
#include <math.h>

#define FREERETURN2 {free_dvector(objtemptemp,0,P);return objtemp;}  /// necessary???
double objfun_mBOO(int i, double *beta, int n,int M, int P,int MCol, double *u, double *dt, int *MNum,double **z, int **MMat, double *xiN, double **B){
    /// for Soh and Huang (2017)
    
    int r,g,gg;
    double lin=0.0,lin2=0.0,objtemp=0.0,elin=0.0, *objtemptemp;
    objtemptemp=dvector(0,P);
    for (r=0;r<=P;r++) objtemptemp[r]=0;
    
    
    ///Get: lin, exp(lin), then objtemp
    
    if (i==0){ /// objective function for i==0
        for (r=0;r<n;r++){
            if (u[i]< dt[r]){
                lin=0.0;
                for (g=0;g<=P;g++) lin+= z[r][g]*beta[g];
                elin=exp(lin);
                objtemp-=xiN[r]*elin;                               /// xiN multiplied
            }
        }
        for (g=0;g<=P;g++){
            for (gg=0;gg<MNum[i];gg++){
                objtemp+=xiN[MMat[i][gg]-1]*z[MMat[i][gg]-1][g]*beta[g];
            }
        }
        
        
        
    }else{ /// objective function for i!=0
        for (r=0;r<n;r++){
            if (u[i]< dt[r]){
                lin=0;
                lin2=0;
                for (g=0;g<=P;g++){
                    lin+= z[r][g]*B[i][g];
                    lin2+= z[r][g]*beta[g];
                }
                elin=exp(lin);
                objtemp-=xiN[r]*exp(lin2);                          /// xiN multiplied
                for (g=0;g<=P;g++){
                    objtemptemp[g]+=xiN[r]*z[r][g]*elin;            /// xiN multiplied
                }
            }
        }
        for (g=0;g<=P;g++){
            objtemp+= objtemptemp[g]*beta[g];
        }
        for (g=0;g<=P;g++){
            for (gg=0;gg<MNum[i];gg++){
                objtemp+=xiN[MMat[i][gg]-1]*z[MMat[i][gg]-1][g]*beta[g];
            }
        }
    }

    /// Objective function is convex. Maximum is where the root exists.
    
    FREERETURN2;
}


