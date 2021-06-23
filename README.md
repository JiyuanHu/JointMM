# JointMM
A joint modeling framework for analyzing zero-inflated longitudinal proportions and time to event data

This package is developed to investigate the association between zero-inflated longitudinal proportions and time to an event,e.g., disease onset. It is in particular applicable to **microbiome proportion data (or called relative abundance data)** from prospective studies. The model JointMM is specifically designed to handle the zero-inflated and highly skewed longitudinal microbial proportion data to examine whether the temporal pattern of microbial presence and/or the non-zero microbial proportions are associated with differences in the time to an event. 

Installation of JointMM in R:

> library("devtools");
> install_github("JiyuanHu/JointMM");

Example code:

> require(JointMM)

> data(dat) #load in the long format of longitudinal and survival data

> zero.prop = mean(dat$Y==0)  

> #1. full model; joint modeling

> res1 = JointMM.func(data=dat,cov.name.long=c('obstime','treatment'),cov.name.surv='treatment',threshold.zero.prop = 0.1,is.longi.model.only = FALSE, quad.n=10) 

> #2. full model; modeling of the longitudinal part only - this is applicable when only longitudinal data is available

> res2 = JointMM.func(data=dat,cov.name.long=c('obstime','treatment'),cov.name.surv='treatment', threshold.zero.prop = 0.1,is.longi.model.only = TRUE, quad.n=10) 

> #3. reduced model-- joint modeling analysis when the proportion of zero RA is lower than the threshold.zero.prop; Here the threshold of 0.6 is just for illustration. For > sample size such like in our demo dataset, a threshold of 10\% is recommended.

> res3 = JointMM.func(data=dat,cov.name.long=c('obstime','treatment'),cov.name.surv='treatment',threshold.zero.prop = 0.6,is.longi.model.only = FALSE, quad.n=10) 

> #4. reduced model-- used when the proportion of zero RA is lower than the threshold.zero.prop; modeling of the longitudinal part only - this is applicable when only longitudinal data is available

> res4 = JointMM.func(data=dat,cov.name.long=c('obstime','treatment'),cov.name.surv='treatment',
                  threshold.zero.prop = 0.6,is.longi.model.only = TRUE, quad.n=10) 

> #names(res1) and names(res3) 

> #c("par.est.Wald","SEs","est.hessian","Wald.Ts","pvals.Wald","status.Wald") 

> #names(res2) and names(res4) 

> #c("par.est.Longonly","SEs","est.hessian","Wald.Ts.Longonly","pvals.Wald.Longonly","status.Wald.Longonly")

> #5. Evaluate the cross-part correlation in the longitudinal sub-model

> res5 = JointMM.func(data=dat,cov.name.long=c('obstime','treatment'),cov.name.surv='treatment',
                    threshold.zero.prop = 0.1,is.longi.model.only = TRUE, cross.part.corr.eval = TRUE,quad.n=10) 


Reference: **Hu J, Wang C, Blaser M, Li H (2021). Joint modeling of zero-inflated longitudinal proportions and time-to-event data with application to a gut microbiome study. Biometrics (Accepted)**
