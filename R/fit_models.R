#' drre: Dynamic Regression with Recurrent Events
#'
#' Recurrent events often arise in follow-up studies where a subject may experience multiple occurrences of the same event. Most regression models with recurrent events tacitly assume constant effects of covariates over time, which may not be realistic in practice. To address time-varying effects, we develop a dynamic regression model to target the mean frequency of recurrent events in Soh and Huang (Biometrics, 2019). Function `drre' fits the regression coefficients in the model based on the proposed estimation procedure which fully exploits observed data. 
#' @param n the number of subjects
#' @param M the number of unique event times in the sample 
#' @param P the number of covariates of interest
#' @param MCol the maximum number of subjects observed in the sample at an event time
#' @param u a vector of ordered unique event times in the sample 
#' @param dt a vector of censoring times
#' @param NNum the number of observations at unique event times; if no tied events, all is 0
#' @param zvec vectorized form of [1 Z]^T matrix. (5 x N)
#' @param Mvec vectorized form of (MCOL x # unique times) ID matrix. ids for unique event times
#' @keywords recurrent events; dynamic regression; mean frequency function;
#' @export
#' @examples
#' drre()

#' @useDynLib drre sandh
drre=function(n,M,P,MCol,u,dt,MNum,zvec,Mvec){
  #Make sure this function gets sorted data in terms of u: Matrix 2  
  temp=.C("sandh",
          as.integer(n),
          as.integer(M),
          as.integer(P),
          as.integer(MCol),
          as.double(u),
          as.double(dt),
          as.integer(MNum), 
          
          as.double(zvec),
          as.integer(Mvec), 
          B_in=as.double(double((M+1)*(P+1))) )
  beta=matrix(temp[["B_in"]],nrow=M+1,ncol=P+1,byrow=T)
  return(list(beta))
}



#' drre_mB: `drre' Based on Nonparametric Multiplier Bootstrap 
#'
#' Function `drre_mB' fits the regression coefficients in the model based on the proposed nonparametric multiplier bootstrap estimation procedure in Soh and Huang (Biometrics, 2019); see related discussion in Rubin (1981), Kosorok (2008) and Huang (2014)
#' @param n the number of subjects
#' @param M the number of unique event times in the sample 
#' @param P the number of covariates of interest
#' @param MCol the maximum number of subjects observed in the sample at an event time
#' @param u a vector of ordered unique event times in the sample 
#' @param dt a vector of censoring times
#' @param NNum the number of observations at unique event times; if no tied events, all is 0
#' @param zvec vectorized form of [1 Z]^T matrix. (5 x N)
#' @param Mvec vectorized form of (MCOL x # unique times) ID matrix. ids for unique event times
#' @param xiN vector of a random sample of size n which following the standard exponential distribution (mean=variance=1)
#' @keywords recurrent events; dynamic regression; mean frequency function;
#' @export
#' @examples
#' drre_mB()

#' @useDynLib drre sandh_mBOO
drre_mB=function(n,M,P,MCol,u,dt,MNum,zvec,Mvec,xiN){
  #MAke sure this function gets sorted data in terms of u
  
  temp=.C("sandh_mBOO",
          as.integer(n),
          as.integer(M),
          as.integer(P),
          as.integer(MCol),
          as.double(u),
          as.double(dt),
          as.integer(MNum),
          
          as.double(zvec),
          as.integer(Mvec),
          as.double(xiN), # N mean 1 var 1 random variables
          B_in=as.double(double((M+1)*(P+1))) )
  beta=matrix(temp[["B_in"]],nrow=M+1,ncol=P+1,byrow=T)
  return(list(beta))
}




#' drre_avg: Average Effects of Covariates Based on `drre' fit
#'
#' Function `drre_ave' provides an average of the estimated effects of covariates coefficients discussed in Soh and Huang (Biometrics, 2019; Section 3.4); see Peng and Huang (1981) for its usage in quantile regression
#' @param t1 Lower limit of the time interval of interest
#' @param t2 Upper limit of the time interval of interest
#' @param result  Data-frame or matrix argument of the vector of unique event times in the sample (increasing order) and the estimated effects of covariates over the event times, i.e., the drre fit
#' @keywords recurrent events; dynamic regression; mean frequency function;
#' @export
#' @examples
#' drre_avg()

drre_avg=function(t1=t1,t2=t2,result=estimated_effects){
  temp.sum=rep(0,dim(result)[2]-1)
  for(i in min(which(result[,1]>t1)):max(which(result[,1]<t2))){
    if (i !=max(which(result[,1]<t2))){
      temp.sum=temp.sum+(result[i+1,1]-result[i,1])*result[i,-1]
    } else {
      temp.sum=temp.sum+(t2-result[i,1])*result[i,-1]
      if (min(which(result[,1]>t1))>1){
        temp.sum=temp.sum+((result[ min(which(result[,1]>t1)),1]-t1)*  #base
                             ((result[min(which(result[,1]>t1)),-1]-result[min(which(result[,1]>t1))-1,-1]) /
                                (result[ min(which(result[,1]>t1)),1]-result[ min(which(result[,1]>t1))-1,1]) *
                                (t1-result[ min(which(result[,1]>t1))-1,1])+result[min(which(result[,1]>t1))-1,-1]) #height
        )
      }                         
    }
  }
  return(temp.sum/(t2-t1))
}


#' drre_TS: A Test Statistic for Evaluating H_0: the effect of a covariate is constant over (t1, t2]
#'
#' Function `drre_TS' provides a test statistic for evaluating the null hypothesis H_0: the effect of a covariate is constant over (t1, t2], which is discussed in Soh and Huang (Biometrics, 2019; Section 3.4).
#' @param t1 Lower limit of the time interval of interest
#' @param t2 Upper limit of the time interval of interest
#' @param result  Data-frame or matrix argument of the vector of unique event times in the sample (increasing order) and the estimated effects of covariates over the event times, i.e., the drre fit
#' @param pwr  The degree of a polynomial weight function
#' @param n the number of subjects in the sample
#' @keywords recurrent events; dynamic regression; mean frequency function;
#' @export
#' @examples
#' drre_avg()

drre_TS=function(t1=t1,t2=t2,result=estimated_effects,pwr=pwr, n=n){
  temp.sum=rep(0,dim(result)[2]-1)  # time+ log(baseline mean) + covariates
  for(i in min(which(result[,1]>t1)):max(which(result[,1]<t2))){
    if (i !=max(which(result[,1]<t2))){
      temp.sum=temp.sum+(result[i+1,1]^(pwr+1)-result[i,1]^(pwr+1))/(pwr+1)*result[i,-1]
    } else {
      temp.sum=temp.sum+(t2^(pwr+1)-result[i,1]^(pwr+1))/(pwr+1)*result[i,-1]
      if (min(which(result[,1]>t1))>1){
        temp.sum=temp.sum+((result[ min(which(result[,1]>t1)),1]^(pwr+1)-t1^(pwr+1))/(pwr+1)*  #base
                             ((result[min(which(result[,1]>t1)),-1]-result[min(which(result[,1]>t1))-1,-1]) /
                                (result[ min(which(result[,1]>t1)),1]-result[ min(which(result[,1]>t1))-1,1]) *
                                (t1-result[ min(which(result[,1]>t1))-1,1])+result[min(which(result[,1]>t1))-1,-1]) #height
        )
      }
    }
  }
  return(temp.sum/(t2-t1)*n^.5)
}
