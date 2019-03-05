#' drre: Dynamic Regression with Recurrent Events
#'
#' Recurrent events often arise in follow-up studies where a subject may experience multiple occurrences of the same event. Most regression models with recurrent events tacitly assume constant effects of covariates over time, which may not be realistic in practice. To address time-varying effects, we develop a dynamic regression model to target the mean frequency of recurrent events in Soh and Huang (Biometrics). Function `drre' fits the regression coefficients in the model based on the proposed estimation procedure, which fully exploits observed data.
#' @param n the number of subjects
#' @param M the number of unique event times in the sample
#' @param P the number of covariates of interest
#' @param MCol the maximum number of subjects observed in the sample at an event time
#' @param u a vector of ordered unique event times in the sample
#' @param dt a vector of censoring times
#' @param MNum the number of observations at unique event times; if no tied events, all is 0
#' @param zvec vectorized form of [1 Z]^T matrix.
#' @param Mvec vectorized form of (MCOL x # unique times) ID matrix. ids for unique event times
#' @keywords recurrent events; dynamic regression; mean frequency function;
#' @export
#' @examples
#' # ------------------
#' # Calling the bladder tumor data (Byar, 1980) from `survival' R package
#' # ------------------
#' require(survival)
#' data("bladder")
#'
#' # ------------------
#' # Data Step
#' # ------------------
#' N=length(unique(bladder1$id))
#' uniq.evnts=sort(unique(bladder1[bladder1$status==1,"stop"]))
#'
#' temp.cen=do.call("rbind",by(bladder1,bladder1$id,tail,n=1))[,c("status","stop")]
#' censoring=ifelse(temp.cen[,"status"]==1,temp.cen[,"stop"]+.1,temp.cen[,"stop"])
#'
#' temp.cov=bladder1[!duplicated(bladder1$id),c("treatment","number","size")]
#' covar=cbind(rep(1,length(unique(bladder1$id))),
#'             as.numeric(temp.cov$treatment=="pyridoxine"),
#'             as.numeric(temp.cov$treatment=="thiotepa"),
#'             temp.cov$number,
#'             temp.cov$size
#' )
#'
#' idMAT=matrix(0,nrow=length(unique(bladder1[bladder1$status==1,"stop"]))
#'              , ncol =  max(table(bladder1[bladder1$status==1,"stop"]))
#' )
#' temp.cnt=1
#' for(i in sort(unique(bladder1[bladder1$status==1,"stop"]))){
#'   temp=sort(bladder1[bladder1$status==1 & bladder1[,"stop"]==i,"id"])
#'   for (j in 1:length(temp)){
#'     idMAT[temp.cnt,j]=temp[j]
#'   }
#'   temp.cnt=temp.cnt+1
#' }
#'
#' # ------------------
#' # Model Fitting
#' # ------------------
#' drre(N,
#'      length(unique(bladder1[bladder1$status==1,"stop"])),
#'      4,
#'      max(table(bladder1[bladder1$status==1,"stop"])),
#'      uniq.evnts,
#'      censoring,
#'      as.vector(table(bladder1[bladder1$status==1,"stop"])),
#'      as.vector(t(covar)),
#'      as.vector(t(idMAT))
#' )
#'

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
#' Function `drre_mB' fits the regression coefficients in the model based on the proposed nonparametric multiplier bootstrap estimation procedure in Soh and Huang (Biometrics); see related discussion in Rubin (1981), Kosorok (2008) and Huang (2014)
#' @param n the number of subjects
#' @param M the number of unique event times in the sample
#' @param P the number of covariates of interest
#' @param MCol the maximum number of subjects observed in the sample at an event time
#' @param u a vector of ordered unique event times in the sample
#' @param dt a vector of censoring times
#' @param MNum the number of observations at unique event times; if no tied events, all is 0
#' @param zvec vectorized form of [1 Z]^T matrix.
#' @param Mvec vectorized form of (MCOL x # unique times) ID matrix. ids for unique event times
#' @param xiN a vector of a random sample of size n from the standard exponential distribution (mean=variance=1)
#' @keywords recurrent events; dynamic regression; mean frequency function;
#' @export
#' @examples
#' # ------------------
#' # Calling the bladder tumor data (Byar, 1980) from `survival' R package
#' # ------------------
#' require(survival)
#' data("bladder")
#'
#' # ------------------
#' # Data Step
#' # ------------------
#' N=length(unique(bladder1$id))
#' uniq.evnts=sort(unique(bladder1[bladder1$status==1,"stop"]))
#'
#' temp.cen=do.call("rbind",by(bladder1,bladder1$id,tail,n=1))[,c("status","stop")]
#' censoring=ifelse(temp.cen[,"status"]==1,temp.cen[,"stop"]+.1,temp.cen[,"stop"])
#'
#' temp.cov=bladder1[!duplicated(bladder1$id),c("treatment","number","size")]
#' covar=cbind(rep(1,length(unique(bladder1$id))),
#'             as.numeric(temp.cov$treatment=="pyridoxine"),
#'             as.numeric(temp.cov$treatment=="thiotepa"),
#'             temp.cov$number,
#'             temp.cov$size
#' )
#'
#' idMAT=matrix(0,nrow=length(unique(bladder1[bladder1$status==1,"stop"]))
#'              , ncol =  max(table(bladder1[bladder1$status==1,"stop"]))
#' )
#' temp.cnt=1
#' for(i in sort(unique(bladder1[bladder1$status==1,"stop"]))){
#'   temp=sort(bladder1[bladder1$status==1 & bladder1[,"stop"]==i,"id"])
#'   for (j in 1:length(temp)){
#'     idMAT[temp.cnt,j]=temp[j]
#'   }
#'   temp.cnt=temp.cnt+1
#' }
#'
#' # ------------------
#' # The Multiplier Bootstrap
#' # ------------------
#' B=1000;set.seed(12345)
#' xiN=matrix(rgamma(N*B,1,1),nrow=N,ncol=B)   # N x B
#' m_results1=list()
#'
#' for (b in 1:B){
#'   m_beta1=drre_mB(N,
#'                   length(unique(bladder1[bladder1$status==1,"stop"])),
#'                   4,
#'                   max(table(bladder1[bladder1$status==1,"stop"])),
#'                   uniq.evnts,
#'                   censoring,
#'                   as.vector(table(bladder1[bladder1$status==1,"stop"])),
#'                   as.vector(t(covar)),
#'                   as.vector(t(idMAT)) ,
#'                   xiN[,b]
#'   )
#'
#'   m_result1=cbind(sort(unique(bladder1[bladder1$status==1,"stop"]))
#'                   ,m_beta1[[1]][-1,])
#'   m_results1[[b]]=m_result1
#' }

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
#' Function `drre_avg' provides an average of the estimated effects of covariates coefficients discussed in Soh and Huang (Biometrics; Section 3.4); see Peng and Huang (2008) for its usage in quantile regression
#' @param t1 Lower limit of the time interval of interest
#' @param t2 Upper limit of the time interval of interest
#' @param result  Data-frame or matrix argument of the vector of unique event times in the sample (increasing order) and the estimated effects of covariates over the event times, i.e., the drre fit
#' @keywords recurrent events; dynamic regression; mean frequency function;
#' @export
#' @examples
#' # ------------------
#' # Calling the bladder tumor data (Byar, 1980) from `survival' R package
#' # ------------------
#' require(survival)
#' data("bladder")
#'
#' # ------------------
#' # Data Step
#' # ------------------
#' N=length(unique(bladder1$id))
#' uniq.evnts=sort(unique(bladder1[bladder1$status==1,"stop"]))
#'
#' temp.cen=do.call("rbind",by(bladder1,bladder1$id,tail,n=1))[,c("status","stop")]
#' censoring=ifelse(temp.cen[,"status"]==1,temp.cen[,"stop"]+.1,temp.cen[,"stop"])
#'
#' temp.cov=bladder1[!duplicated(bladder1$id),c("treatment","number","size")]
#' covar=cbind(rep(1,length(unique(bladder1$id))),
#'             as.numeric(temp.cov$treatment=="pyridoxine"),
#'             as.numeric(temp.cov$treatment=="thiotepa"),
#'             temp.cov$number,
#'             temp.cov$size
#' )
#'
#' idMAT=matrix(0,nrow=length(unique(bladder1[bladder1$status==1,"stop"]))
#'              , ncol =  max(table(bladder1[bladder1$status==1,"stop"]))
#' )
#' temp.cnt=1
#' for(i in sort(unique(bladder1[bladder1$status==1,"stop"]))){
#'   temp=sort(bladder1[bladder1$status==1 & bladder1[,"stop"]==i,"id"])
#'   for (j in 1:length(temp)){
#'     idMAT[temp.cnt,j]=temp[j]
#'   }
#'   temp.cnt=temp.cnt+1
#' }
#'
#' # ------------------
#' # Model Fitting
#' # ------------------
#' BETA1=
#'      drre(N,
#'      length(unique(bladder1[bladder1$status==1,"stop"])),
#'      4,
#'      max(table(bladder1[bladder1$status==1,"stop"])),
#'      uniq.evnts,
#'      censoring,
#'      as.vector(table(bladder1[bladder1$status==1,"stop"])),
#'      as.vector(t(covar)),
#'      as.vector(t(idMAT))
#' )
#' trim.result=cbind(uniq.evnts,BETA1[[1]][-1,])
#'
#' # ------------------
#' # Estimates of Average Effects of Covariates in time interval (t1,t2]
#' # ------------------
#' t1=5; t2=53
#' drre_avg(t1,t2,trim.result)[-1]
#'

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
#' Function `drre_TS' provides a test statistic for evaluating the null hypothesis H_0: the effect of a covariate is constant over (t1, t2], which is discussed in Soh and Huang (Biometrics; Section 3.4).
#' @param t1 Lower limit of the time interval of interest
#' @param t2 Upper limit of the time interval of interest
#' @param result  Data-frame or matrix argument of the vector of unique event times in the sample (increasing order) and the estimated effects of covariates over the event times, i.e., the drre fit
#' @param pwr  The degree of a polynomial weight function
#' @param n The number of subjects in the sample
#' @keywords recurrent events; dynamic regression; mean frequency function;
#' @export
#' @examples
#' # ------------------
#' # Calling the bladder tumor data (Byar, 1980) from `survival' R package
#' # ------------------
#' require(survival)
#' data("bladder")
#'
#' # ------------------
#' # Data Step
#' # ------------------
#' N=length(unique(bladder1$id))
#' uniq.evnts=sort(unique(bladder1[bladder1$status==1,"stop"]))
#'
#' temp.cen=do.call("rbind",by(bladder1,bladder1$id,tail,n=1))[,c("status","stop")]
#' censoring=ifelse(temp.cen[,"status"]==1,temp.cen[,"stop"]+.1,temp.cen[,"stop"])
#'
#' temp.cov=bladder1[!duplicated(bladder1$id),c("treatment","number","size")]
#' covar=cbind(rep(1,length(unique(bladder1$id))),
#'             as.numeric(temp.cov$treatment=="pyridoxine"),
#'             as.numeric(temp.cov$treatment=="thiotepa"),
#'             temp.cov$number,
#'             temp.cov$size
#' )
#'
#' idMAT=matrix(0,nrow=length(unique(bladder1[bladder1$status==1,"stop"]))
#'              , ncol =  max(table(bladder1[bladder1$status==1,"stop"]))
#' )
#' temp.cnt=1
#' for(i in sort(unique(bladder1[bladder1$status==1,"stop"]))){
#'   temp=sort(bladder1[bladder1$status==1 & bladder1[,"stop"]==i,"id"])
#'   for (j in 1:length(temp)){
#'     idMAT[temp.cnt,j]=temp[j]
#'   }
#'   temp.cnt=temp.cnt+1
#' }
#'
#' # ------------------
#' # Model Fitting
#' # ------------------
#' BETA1=
#'      drre(N,
#'      length(unique(bladder1[bladder1$status==1,"stop"])),
#'      4,
#'      max(table(bladder1[bladder1$status==1,"stop"])),
#'      uniq.evnts,
#'      censoring,
#'      as.vector(table(bladder1[bladder1$status==1,"stop"])),
#'      as.vector(t(covar)),
#'      as.vector(t(idMAT))
#' )
#' trim.result=cbind(uniq.evnts,BETA1[[1]][-1,])
#'
#' #---------------------------
#' # Test Statistic for 'Constant' Effect
#' #---------------------------
#' t1=5; t2=53
#' pwr=1
#'
#' trim.result_ConstantEffect=cbind(uniq.evnts,t(t(trim.result[,-1])-drre_avg(t1,t2,trim.result)))
#' drre_TS(t1,t2,trim.result_ConstantEffect,pwr=pwr,n=N)
#'

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
