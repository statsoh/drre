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
