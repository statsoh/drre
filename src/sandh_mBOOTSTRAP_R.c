#include <R.h>
#include "nrutilR.h"
#include "NR_mBOOTSTRAP_R.h"

void sandh_mBOO(int *n_in, int *M_in, int *P_in,int *MCol_in, double *u_in, double *d_in,int *MNum_in, double *zvec_in , int *MMat_in ,double *xiN_in, double *B_in){
    int i,j,g;

    int n= n_in[0];
    int M= M_in[0];
    int P= P_in[0];
    int MCol= MCol_in[0];
    
    int *MNum=MNum_in;
    int  **MMat;
    double *u=u_in,*dt=d_in, *xiN=xiN_in ,**z, **B;
    
    z=dmatrix(0,n-1,0,P);
    MMat=imatrix(0,M-1,0,MCol-1);
    B=dmatrix(0,M,0,P);
    for(i=0;i<n;i++) z[i]= &zvec_in[i*(P+1)];
    for(i=0;i<M;i++) MMat[i]= &MMat_in[i*MCol];
    for(i=0;i<M+1;i++) B[i]= &B_in[i*(P+1)];
#if 0
    /// ///////////////////////////////////////////////////////DATA checking in R :Start
    Rprintf("Data checking in R: Start\n");
    Rprintf("#total=%d\t#unique-event-time=%d\t#covariates=%d\t#maximum-possible at m\n",n,M,P,MCol);

    for(i=0;i<M;i++) Rprintf("%lf ",u[i]);Rprintf("\n\n");
    for(i=0;i<M;i++) Rprintf("%d ",MNum[i]);Rprintf("\n\n");
    for(i=0;i<n;i++) Rprintf("%lf ",dt[i]);Rprintf("\n\n"); /// last T in each observation
    for(i=0;i<n;i++) Rprintf("%lf ",xiN[i]);Rprintf("\n\n");

    for (i=0;i<n;i++){
        for (j=0;j<=P;j++){
            Rprintf("%lf ",z[i][j]);
        }Rprintf("\n");
    }Rprintf("\n");
    for (i=0;i<M;i++){
        for (j=0;j<MCol;j++){
            Rprintf("%lf ",MMat[i][j]);
        }Rprintf("\n");
    }Rprintf("\n");
    ///    for (i=0;i<M+1;i++){
    ///        for (j=0;j<=P;j++){
    ///            Rprintf("%lf ",B[i][j]);
    ///        }Rprintf("\n");
    ///    }
    Rprintf("Data checking in C: done\n");
    /// ///////////////////////////////////////////////////////DATA checking in C : Done
#endif
    /// ///////////////////////////////////////////////////////Estimation of Beta :Start
    
    double *tempbeta=(double *)malloc((P+1)*sizeof(double));
    
    for (i=0;i<M;i++){
        ///       Rprintf("\n================================================================================");
        ///        Rprintf("\n%d-th unique event time, out of %d;\t Total # of ind. %d\n",i+1,M,n);
        
        ///        Rprintf("------------Estimates of beta from the previous step; for i-th EE------------\n");
        ///        for (g=0;g<=P;g++) Rprintf("%lf\t",B[i][g]);
        ///        Rprintf("\n----------------------------------------------------------------\n");
        
        tempbeta=NR_mBOO(i,n,M,P,MCol,u,dt,MNum,z,MMat,xiN,B);
        for (g=0;g<=P;g++) B[i+1][g]=tempbeta[g];
        
    }
    
    /// Print final Beta[]
    for (i=0;i<=M;i++){
        for (j=0;j<=P;j++){
            ///            Rprintf("\t%lf",B[i][j]);
            B_in[i*(P+1)+j]=B[i][j];
        }//Rprintf("\n");
    }
    
    
}

