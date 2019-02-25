#include <R.h>
#include "nrutilR.h"
#include <math.h>
#include "ludecompR.h"
#include "userfunR.h"
#include "objfunR.h"

#define NTRIAL (int)100
#define TOLX 1.0e-8
#define TOLF 1.0e-10
#define TOLdB 1.5 /// 150

#define FREERETURN {free_dmatrix(fjac,1,P+1,1,P+1);free_dvector(fvec,1,P+1);free_dvector(tempdbeta,1,P+1);free_ivector(indx,1,P+1);free_dvector(tempBeta,0,P);return beta;}

double *NR(int i, int n, int M, int P,int MCol, double *u, double *dt, int *MNum ,double **z, int **MMat,double **B)
{ /// sandh 2017

    static int td=1;
    static double *beta=NULL; /// 'static' is important to ensure this NR return the correct address of beta.
    if (td) {beta=dvector(0,P);
            td=0;}
    double errf,errx,d,*tempdbeta,*fvec,**fjac,big,temp,obj_f1,obj_f2,*tempBeta;
    int k,g,gg,*indx;

    tempBeta=dvector(0,P);

    indx=ivector(1,P+1);        /// 'Numerical Recipes' functions --- memory allocation
    tempdbeta=dvector(1,P+1);
    fvec=dvector(1,P+1);
    fjac=dmatrix(1,P+1,1,P+1);

    for (g=0;g<=P;g++) beta[g]=B[i][g];
    for (k=1;k<=NTRIAL;k++){
        /// In userfun, fvec and fjac are initialized and calculated at given beta: beta={0,0,0,0,0} when k=1.
        userfun(i,beta,fvec,fjac,n,M,P,MCol,u,dt,MNum,z,MMat,B); /// check done
#if 1
        Rprintf("\n%d-th unique obs event(EE); %d-th NR iter:",i+1,k);
        Rprintf("\n------------beta (input) -------------\n");
        for (g=0;g<=P;g++) Rprintf("%.10lf\t",beta[g]);
        Rprintf("\n------------fvec (output;input beta)-------------\n");
        for (g=0;g<=P;g++) Rprintf("%.10lf\t",fvec[g+1]);
        Rprintf("\n------------fjac (output;input beta)-------------\n");
        for (g=1;g<=P+1;g++){
            for (gg=1;gg<=P+1;gg++){
                Rprintf("%.10lf\t",fjac[g][gg]);
            }Rprintf("\n");
        }
#endif
        for (g=1;g<=P+1;g++){
            big=0.0;
            for (gg=1;gg<=P+1;gg++){
                big+=fabs(fjac[gg][g]); /// Check whether a column has all zeros.
                if (big==0.0) {
                    warning("Singular jacobian matrix (fjac) in routine N-R\n");
                    FREERETURN; /// return beta[] now
                }
            }
        }

        errf=0.0;
        for (g=0;g<=P;g++) errf+=fabs(fvec[g+1]);
        if (errf==0.0) {
            Rprintf("\nManhattan distance of Estimating Equation is close enough to zero, which is good. Go to next EE.\n");
            FREERETURN;
        }
     Rprintf("\nManhattan distance of EE : \t%.4lf",errf);
        for (g=0;g<=P;g++) tempdbeta[g+1]= -fvec[g+1];
        ludcmp(fjac,P+1,indx,&d);   /// LU decomposition: fjac changes
        lubksb(fjac,P+1,indx,tempdbeta); ///Backward substitution using changed fjac: (fjac)*(dbeta)=(-fvec) /// answer (dbeta) is given back to (tempdbeta)

        errx=0.0;
        for(g=0;g<=P;g++){
            errx+= fabs(tempdbeta[g+1]);
            tempBeta[g] = beta[g];
            tempBeta[g] += tempdbeta[g+1];
        }

        if (errx==0.0) {
                beta[g] += tempdbeta[g+1];
                Rprintf("\nSmall enough changes in beta, which is good.\n");
                FREERETURN;  /// small enough changes in beta? Then returns beta[], which is good.
        }
        Rprintf("\nDistance(dBeta;errx_before_WHILE)=%.8lf",errx);

/// The Concave Objective function

        obj_f1=objfun(i,beta,n,M,P,MCol,u,dt,MNum,z,MMat,B);       ///initial values at the N-R step
        obj_f2=objfun(i,tempBeta,n,M,P,MCol,u,dt,MNum,z,MMat,B);    ///updated values : objfun should be greater than obj_f1.
        userfun(i,tempBeta,fvec,fjac,n,M,P,MCol,u,dt,MNum,z,MMat,B);
        errf=0.0;
        for (g=0;g<=P;g++) errf+=fabs(fvec[g+1]);
        
        gg=0;
        while( ( (obj_f1 > obj_f2) && (gg < NTRIAL ) )  ||  ( (errx >=TOLdB*(P+1)) && (gg < NTRIAL ) )  ){
            temp=0;  gg++;
            for(g=0;g<=P;g++){
                tempBeta[g] = beta[g]+tempdbeta[g+1]/pow(2,gg);  ///
                temp+= fabs(tempdbeta[g+1]/pow(2,gg));         /// for some reasons, 2^gg gives error in LUdecomp... pow(2,gg).;;;;;
            }
            errx=temp;
            obj_f2=objfun(i,tempBeta,n,M,P,MCol,u,dt,MNum,z,MMat,B);
            Rprintf("\n%d-th halving; \t%lf (increment in OBJ);\t then Distance(dBeta)=%lf",gg,obj_f2-obj_f1,errx);
        }
        Rprintf("\n%lf (obj_f1);\t %lf (obj_f2) ----- difference: %lf\nDistance(dBeta)=%lf \t %d-th halving\n",obj_f1,obj_f2,obj_f2-obj_f1,errx,gg);
        
    ///////// NEW additional criteria
        userfun(i,tempBeta,fvec,fjac,n,M,P,MCol,u,dt,MNum,z,MMat,B);
        errf=0.0;
        for (g=0;g<=P;g++) errf+=fabs(fvec[g+1]);
        if (errf<TOLF) {
            Rprintf("\nManhattan distance of Estimating Equation is close enough to zero, which is good. Go to next EE.\n");
            for (g=0;g<=P;g++) beta[g]=tempBeta[g];
            FREERETURN;
        }
    /////////////////////////////////////
        if (obj_f1 > obj_f2){
            Rprintf("\nWarning!! No improvement on Objective function value. Return the initial beta[] (not want to see); \n# iter: %d ",gg);
            warning("\nWarning!! No improvement on Objective function value. Return the initial beta[] (not want to see); \n# iter: %d ",gg);
            FREERETURN;
        }else
        {
                        Rprintf("%lf (increment in OBJ); Adding half of dBeta for %d times ; %d-th NR iteration\n",obj_f2-obj_f1,gg,k);
            for (g=0;g<=P;g++) beta[g]=tempBeta[g];
            if (fabs(obj_f2-obj_f1)< TOLF) {
                Rprintf("\n-------------------------------Distance(dBeta)=%.10lf \t M.Dis(Beta)=%.10lf\n(OBJ2-OBJ1)>=0 is close enough to zero, which is good. Go to next EE. \n",errx,errf);
                FREERETURN;}
            continue;
        }
    }
    FREERETURN;

}

