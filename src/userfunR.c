#include <R.h>
#include <math.h>
#include "nrutilR.h"


/// sandh
void userfun(int i, double *beta, double *fvec, double **fjac, int n,int M, int P,int MCol, double *u ,double *dt, int *MNum,double **z,int **MMat, double **B){ ///Form *fvec and **fjac at beta[]

    int g,r,gg;
    double lin=0.0,lin2=0.0,elin=0.0,elin2=0.0;

    //fvec and fjac initialization to zero
    for (g=0;g<=P;g++){
        fvec[g+1]=0;
        for (gg=0;gg<=P;gg++) fjac[g+1][gg+1]=0;
    }

    //Get: fvec and fjac
    if (i==0){
        
        for (r=0;r<n;r++){
            if (u[i] < dt[r] ){ /// uses not censored r... /// u : unique time m. (M x 1)  ; dt : Censoring time (n x 1)
                lin=0.0;
                for (g=0;g<=P;g++) lin+= z[r][g]*beta[g];
                elin=exp(lin);
                for (g=0;g<=P;g++){
                    fvec[g+1]-=z[r][g]*elin;
                    for (gg=0;gg<=P;gg++){
                        fjac[g+1][gg+1]-=z[r][g]*z[r][gg]*elin;
                    }
                }
            }
        }
        for (g=0;g<=P;g++){
            for (gg=0;gg<MNum[i];gg++){
                fvec[g+1]+=z[MMat[i][gg]-1][g];  /// MMat[i][gg] means (gg+1)-th #id among the covariates whose event time is u[i]
            }
        }
    }else{
        for (r=0;r<n;r++){
            if (u[i] < dt[r] ){ /// uses not censored r
                lin=0.0;
                lin2=0.0;
                for (g=0;g<=P;g++){
                    lin+= z[r][g]*beta[g];
                    lin2+= z[r][g]*B[i][g]; /// We will save beta estimate to B[i+1][g]
                }
                elin=exp(lin);
                elin2=exp(lin2);
                for (g=0;g<=P;g++){
                    fvec[g+1] -= z[r][g]*(elin-elin2);
                    for (gg=0;gg<=P;gg++){
                        fjac[g+1][gg+1]-=z[r][g]*z[r][gg]*elin;
                    }
                }
            }
        }
        for (g=0;g<=P;g++){
            for (gg=0;gg<MNum[i];gg++){
                fvec[g+1]+=z[MMat[i][gg]-1][g];
            }
        }//// now we have fjac and fvec at global memory
    }
}
