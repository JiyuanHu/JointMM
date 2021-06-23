\name{JointMM-package}
\alias{JointMM-package}
\alias{JointMM}
\docType{package}
\title{
\packageTitle{JointMM}
}
\description{
\packageDescription{JointMM}
}
\details{

A joint modeling tool for analyzing zero-inflated longitudinal proportions and time-to-event data
\packageDESCRIPTION{JointMM}
\packageIndices{JointMM}
This package is developed to investigate the association between zero-inflated longitudinal proportions and time to an event,e.g., disease onset. It is in particular applicable to microbiome proportion data from prospective studies,  specifically designed to handle the zero-inflated and highly skewed longitudinal microbial proportion data to examine whether the temporal pattern of microbial presence and/or the non-zero microbial proportions are associated with differences in the time to an event. 
The demo dataset `dat` illustrates how to organize the longitudinal and survival data for the joint modeling analysis. The main function of this package `JointMM.func()`. In the example we illustrated how to fit JointMM to the demo data to perform joint modeling analysis.
}
\author{
\packageAuthor{JointMM}
Jiyuan Hu
Maintainer: Jiyuan Hu\packageMaintainer{JointMM}
}
\references{
Hu J, Wang C, Blaser M, Li H (2021). Joint modeling of zero-inflated longitudinal proportions and time-to-event data with application to a gut microbiome study. Biometrics(Accepted).
}
\keyword{ joint modeling, microbiome, zero inflation, skewness }
\examples{
require(JointMM)
data(dat) #load in the long format of longitudinal and survival data

# the proportion of zero RAs for the demo dataset is 52%. We recommend that a reduced model of JointMM is fitted when the zero proportion is too low (for example, <10%) and this threshold depends on the sample size and the number of observations for the data
zero.prop = mean(dat$Y==0)  

#1. full model; joint modeling
res1 = JointMM.func(data=dat,cov.name.long=c('obstime','treatment'),cov.name.surv='treatment',
                  threshold.zero.prop = 0.1,is.longi.model.only = FALSE, quad.n=10) 
#2. full model; modeling of the longitudinal part only - this is applicable when only longitudinal data is available
res2 = JointMM.func(data=dat,cov.name.long=c('obstime','treatment'),cov.name.surv='treatment',
                  threshold.zero.prop = 0.1,is.longi.model.only = TRUE, quad.n=10) 
#3. reduced model-- joint modeling analysis when the proportion of zero RA is lower than the threshold.zero.prop; Here the threshold of 0.6 is just for illustration. For sample size such like in our demo dataset, a threshold of 10\% is recommended.
res3 = JointMM.func(data=dat,cov.name.long=c('obstime','treatment'),cov.name.surv='treatment',
                  threshold.zero.prop = 0.6,is.longi.model.only = FALSE, quad.n=10) 

#4. reduced model-- used when the proportion of zero RA is lower than the threshold.zero.prop; modeling of the longitudinal part only - this is applicable when only longitudinal data is available
res4 = JointMM.func(data=dat,cov.name.long=c('obstime','treatment'),cov.name.surv='treatment',
                  threshold.zero.prop = 0.6,is.longi.model.only = TRUE, quad.n=10) 
#names(res1) and names(res3) 
#c("par.est.Wald","SEs","est.hessian","Wald.Ts","pvals.Wald","status.Wald") 

#names(res2) and names(res4) 
#c("par.est.Longonly","SEs","est.hessian","Wald.Ts.Longonly","pvals.Wald.Longonly","status.Wald.Longonly")

#5. Evaluate the cross-part correlation in the longitudinal sub-model
res5 = JointMM.func(data=dat,cov.name.long=c('obstime','treatment'),cov.name.surv='treatment',
                    threshold.zero.prop = 0.1,is.longi.model.only = TRUE, cross.part.corr.eval = TRUE,quad.n=10) 

}
