# drre
This package provides useful R functions for the proposed estimation and inference procedures, discussed in `Dynamic Regression with Recurrent Events' of Jay Soh and Dr. Eugene Huang.

Email: statsoh@gmail.com (Jay Soh)


# To install this package
```
library(devtools)
install_github("statsoh/drre")
```

# drre::drre()
Function `drre' fits time-varying regression coefficients of the proposed dynamic regression model with recurrent events. The  estimation procedure fully exploits observed data.
### Usage
```
drre(n, M, P, MCol, u, dt, MNum, zvec, Mvec)
```
### Arguments
- n:     the number of subjects
- M:     the number of unique event times in the sample
- P:     the number of covariates of interest
- MCol:     the maximum number of subjects observed in the sample at an event time
- u:     a vector of ordered unique event times in the sample
- dt:     a vector of censoring times
- MNum:     the number of observations at unique event times; if no tied events, all is 0
- zvec:     vectorized form of [1 Z]^T matrix.
- Mvec:     vectorized form of (MCOL x # unique times) ID matrix. ids for unique event times

### Example
```
# ------------------
# Calling the bladder tumor data (Byar, 1980) from `survival' R package
# ------------------
require(survival)
data("bladder")

# ------------------
# Data Step
# ------------------
N=length(unique(bladder1$id))
uniq.evnts=sort(unique(bladder1[bladder1$status==1,"stop"]))

temp.cen=do.call("rbind",by(bladder1,bladder1$id,tail,n=1))[,c("status","stop")]
censoring=ifelse(temp.cen[,"status"]==1,temp.cen[,"stop"]+.1,temp.cen[,"stop"])

temp.cov=bladder1[!duplicated(bladder1$id),c("treatment","number","size")]
covar=cbind(rep(1,length(unique(bladder1$id))),
            as.numeric(temp.cov$treatment=="pyridoxine"),
            as.numeric(temp.cov$treatment=="thiotepa"),
            temp.cov$number,
            temp.cov$size
)

idMAT=matrix(0,nrow=length(unique(bladder1[bladder1$status==1,"stop"]))
             , ncol =  max(table(bladder1[bladder1$status==1,"stop"]))
)
temp.cnt=1
for(i in sort(unique(bladder1[bladder1$status==1,"stop"]))){
  temp=sort(bladder1[bladder1$status==1 & bladder1[,"stop"]==i,"id"])
  for (j in 1:length(temp)){
    idMAT[temp.cnt,j]=temp[j]
  }
  temp.cnt=temp.cnt+1
}

# ------------------
# Model Fitting
# ------------------
drre(N,
     length(unique(bladder1[bladder1$status==1,"stop"])),
     4,
     max(table(bladder1[bladder1$status==1,"stop"])),
     uniq.evnts,
     censoring,
     as.vector(table(bladder1[bladder1$status==1,"stop"])),
     as.vector(t(covar)),
     as.vector(t(idMAT))
)
```

# drre::drre_mB()
Function ‘drre_mB’ fits the regression coefficients in the model based on the proposed nonparametric multiplier bootstrap estimation procedure in Soh and Huang (Biometrics); see related discussion in Rubin (1981), Kosorok (2008) and Huang (2014)

### Usage
```
drre_mB(n, M, P, MCol, u, dt, MNum, zvec, Mvec, xiN)
```
### Arguments
- n:     the number of subjects
- M:     the number of unique event times in the sample
- P:     the number of covariates of interest
- MCol:     the maximum number of subjects observed in the sample at an event time
- u:     a vector of ordered unique event times in the sample
- dt:     a vector of censoring times
- MNum:     the number of observations at unique event times; if no tied events, all is 0
- zvec:     vectorized form of [1 Z]^T matrix.
- Mvec:     vectorized form of (MCOL x # unique times) ID matrix. ids for unique event times
- xiN:	a vector of a random sample of size n from the standard exponential distribution (mean=variance=1)

### Example
```
# ------------------
# Calling the bladder tumor data (Byar, 1980) from `survival' R package
# ------------------
require(survival)
data("bladder")

# ------------------
# Data Step
# ------------------
N=length(unique(bladder1$id))
uniq.evnts=sort(unique(bladder1[bladder1$status==1,"stop"]))

temp.cen=do.call("rbind",by(bladder1,bladder1$id,tail,n=1))[,c("status","stop")]
censoring=ifelse(temp.cen[,"status"]==1,temp.cen[,"stop"]+.1,temp.cen[,"stop"])

temp.cov=bladder1[!duplicated(bladder1$id),c("treatment","number","size")]
covar=cbind(rep(1,length(unique(bladder1$id))),
            as.numeric(temp.cov$treatment=="pyridoxine"),
            as.numeric(temp.cov$treatment=="thiotepa"),
            temp.cov$number,
            temp.cov$size
)

idMAT=matrix(0,nrow=length(unique(bladder1[bladder1$status==1,"stop"]))
             , ncol =  max(table(bladder1[bladder1$status==1,"stop"]))
)
temp.cnt=1
for(i in sort(unique(bladder1[bladder1$status==1,"stop"]))){
  temp=sort(bladder1[bladder1$status==1 & bladder1[,"stop"]==i,"id"])
  for (j in 1:length(temp)){
    idMAT[temp.cnt,j]=temp[j]
  }
  temp.cnt=temp.cnt+1
}

# ------------------
# The Multiplier Bootstrap
# ------------------
B=1000;set.seed(12345)
xiN=matrix(rgamma(N*B,1,1),nrow=N,ncol=B)   # N x B
m_results1=list()

for (b in 1:B){
  m_beta1=drre_mB(N,
                  length(unique(bladder1[bladder1$status==1,"stop"])),
                  4,
                  max(table(bladder1[bladder1$status==1,"stop"])),
                  uniq.evnts,
                  censoring,
                  as.vector(table(bladder1[bladder1$status==1,"stop"])),
                  as.vector(t(covar)),
                  as.vector(t(idMAT)) ,
                  xiN[,b]
  )

  m_result1=cbind(sort(unique(bladder1[bladder1$status==1,"stop"]))
                  ,m_beta1[[1]][-1,])
  m_results1[[b]]=m_result1
}
```

# drre::drre_avg()
Function ‘drre_avg’ provides an average of the estimated effects of covariates coefficients discussed in Soh and Huang (Biometrics; Section 3.4); see Peng and Huang (2008) for its usage in quantile regression
### Usage
```
drre_avg(t1, t2, result)
```
### Arguments
- t1:       Lower limit of the time interval of interest
- t2:       Upper limit of the time interval of interest
- result:   Data-frame or matrix argument of the vector of unique event times in the sample (increasing order) and the estimated effects of covariates over the event times, i.e., `drre' fit

### Example
```
# ------------------
# Calling the bladder tumor data (Byar, 1980) from `survival' R package
# ------------------
require(survival)
data("bladder")

# ------------------
# Data Step
# ------------------
N=length(unique(bladder1$id))
uniq.evnts=sort(unique(bladder1[bladder1$status==1,"stop"]))

temp.cen=do.call("rbind",by(bladder1,bladder1$id,tail,n=1))[,c("status","stop")]
censoring=ifelse(temp.cen[,"status"]==1,temp.cen[,"stop"]+.1,temp.cen[,"stop"])

temp.cov=bladder1[!duplicated(bladder1$id),c("treatment","number","size")]
covar=cbind(rep(1,length(unique(bladder1$id))),
            as.numeric(temp.cov$treatment=="pyridoxine"),
            as.numeric(temp.cov$treatment=="thiotepa"),
            temp.cov$number,
            temp.cov$size
)

idMAT=matrix(0,nrow=length(unique(bladder1[bladder1$status==1,"stop"]))
             , ncol =  max(table(bladder1[bladder1$status==1,"stop"]))
)
temp.cnt=1
for(i in sort(unique(bladder1[bladder1$status==1,"stop"]))){
  temp=sort(bladder1[bladder1$status==1 & bladder1[,"stop"]==i,"id"])
  for (j in 1:length(temp)){
    idMAT[temp.cnt,j]=temp[j]
  }
  temp.cnt=temp.cnt+1
}

# ------------------
# Model Fitting
# ------------------
BETA1=
     drre(N,
     length(unique(bladder1[bladder1$status==1,"stop"])),
     4,
     max(table(bladder1[bladder1$status==1,"stop"])),
     uniq.evnts,
     censoring,
     as.vector(table(bladder1[bladder1$status==1,"stop"])),
     as.vector(t(covar)),
     as.vector(t(idMAT))
)
trim.result=cbind(uniq.evnts,BETA1[[1]][-1,])

# ------------------
# Estimates of Average Effects of Covariates in time interval (t1,t2]
# ------------------
t1=5; t2=53
drre_avg(t1,t2,trim.result)[-1]
```

# drre::drre_TS()
Function ‘drre_TS’ provides a test statistic for evaluating the null hypothesis H_0: the effect of a covariate is constant over (t1, t2], which is discussed in Soh and Huang (Biometrics; Section 3.4).
### Usage
```
drre_TS(t1, t2, result, pwr, n)
```
### Arguments
- t1:       Lower limit of the time interval of interest
- t2:       Upper limit of the time interval of interest
- result:   Data-frame or matrix argument of the vector of unique event times in the sample (increasing order) and the estimated effects of covariates over the event times, i.e., `drre' fit
- pwr:       The degree of a polynomial weight function
- n:        the number of subjects in the sample
### Example
```
# ------------------
# Calling the bladder tumor data (Byar, 1980) from `survival' R package
# ------------------
require(survival)
data("bladder")

# ------------------
# Data Step
# ------------------
N=length(unique(bladder1$id))
uniq.evnts=sort(unique(bladder1[bladder1$status==1,"stop"]))

temp.cen=do.call("rbind",by(bladder1,bladder1$id,tail,n=1))[,c("status","stop")]
censoring=ifelse(temp.cen[,"status"]==1,temp.cen[,"stop"]+.1,temp.cen[,"stop"])

temp.cov=bladder1[!duplicated(bladder1$id),c("treatment","number","size")]
covar=cbind(rep(1,length(unique(bladder1$id))),
            as.numeric(temp.cov$treatment=="pyridoxine"),
            as.numeric(temp.cov$treatment=="thiotepa"),
            temp.cov$number,
            temp.cov$size
)

idMAT=matrix(0,nrow=length(unique(bladder1[bladder1$status==1,"stop"]))
             , ncol =  max(table(bladder1[bladder1$status==1,"stop"]))
)
temp.cnt=1
for(i in sort(unique(bladder1[bladder1$status==1,"stop"]))){
  temp=sort(bladder1[bladder1$status==1 & bladder1[,"stop"]==i,"id"])
  for (j in 1:length(temp)){
    idMAT[temp.cnt,j]=temp[j]
  }
  temp.cnt=temp.cnt+1
}

# ------------------
# Model Fitting
# ------------------
BETA1=
     drre(N,
     length(unique(bladder1[bladder1$status==1,"stop"])),
     4,
     max(table(bladder1[bladder1$status==1,"stop"])),
     uniq.evnts,
     censoring,
     as.vector(table(bladder1[bladder1$status==1,"stop"])),
     as.vector(t(covar)),
     as.vector(t(idMAT))
)
trim.result=cbind(uniq.evnts,BETA1[[1]][-1,])

#---------------------------
# Test Statistic for 'Constant' Effect
#---------------------------
t1=5; t2=53
pwr=1

trim.result_ConstantEffect=cbind(uniq.evnts,t(t(trim.result[,-1])-drre_avg(t1,t2,trim.result)))
drre_TS(t1,t2,trim.result_ConstantEffect,pwr=pwr,n=N)

```


