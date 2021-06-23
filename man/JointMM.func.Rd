\name{JointMM.func}
\alias{JointMM.func}
\title{
The main function of the joint modeling tool.
}
\description{This function allows investigators to jointly analyze zero-inflated longitudinal proportion data (typically microbiome relative abundance data) and survival data to examine whether the temporal pattern of microbial presence and/or the non-zero microbial proportions are associated with differences in the time to an event. When only longitudinal data is available, this function also allows the investigation of the association between (time-varying) covariates and longitudinal proportions.
}
\usage{
#joint modeling of longitudinal and survival data
JointMM.func(data, cov.name.long, cov.name.surv, threshold.zero.prop = 0.1, is.longi.model.only = FALSE, cross.part.corr.eval = FALSE, quad.n = 10)
#modeling of longitudinal data only
JointMM.func(data, cov.name.long, cov.name.surv, threshold.zero.prop = 0.1, is.longi.model.only = TRUE, cross.part.corr.eval = FALSE,quad.n = 10)
}

\arguments{
  \item{data}{
A data frame containing the longitudinal microbial proportions of one taxon (or other types of zero-inflated proportion data), the survival outcome, and (time-dependent) covariates for each subject in the long format. The demo dataset of this package `dat` illustrates how to organize the longitudinal and survival data for the joint modeling analysis
}
  \item{cov.name.long}{
A vector of covariates to be adjusted in the longitudinal sub-model. `cov.name.long` should be consistent with that of the colnames of data
}
  \item{cov.name.surv}{
A vector of covariates to be adjusted in the survival sub-model. `cov.name.surv` should be consistent with that of the colnames of data
}
  \item{threshold.zero.prop}{
The threshold of zero proportions to indicate whether to fit the reduced model of JointMM. A full model is fitted if the proportion of zero RA of the data is higher than `threshold.zero.prop`, and a reduced model if fitted otherwise. The threshold depends on the sample size and number of observations of the data. 
}
  \item{is.longi.model.only}{
An indicator of whether only the longitudinal sub-model of JointMM is used. `is.longi.model.only`= TRUE if only longitudinal data is available
}
  \item{cross.part.corr.eval}{
An indicator of whether the cross-part correlation is evaluated for the longitudinal sub-model of JointMM. We do not recommend assessing the cross-part correlation for conducting microbiome-wide association screening and therefore the default value is set as FALSE. 
}
  \item{quad.n}{
The number of points used for numerical integration under the Gauss-Hermite rule. `quad.n` = 10 is set by default in the package
}
}

\value{
  \item{par.est.Wald}{A vector containing the point estimates of parameters in the joint model or in the longitudinal sub-model when `is.longi.model.only` = TRUE}
  \item{SEs}{A vector containing the estimated standard error of the point estimates in the joint model or in the longitudinal sub-model when `is.longi.model.only` = TRUE}
  \item{est.hessian}{The estimated hessian matrix}
  \item{Wald.Ts}{A vector containing the Wald statistics for testing the global null (Wald.all), the presence part (Wald.delta2), and the non-zero proportion part (Wald.delta2), respectively for joint modeling analysis. Only Wald.delta2 is available if fitting reduced model of JointMM}
  \item{pvals.Wald}{The corresponding p-values of the Wald statistics Wald.Ts for joint modeling analysis}
  \item{status.Wald}{An indicator of whether the optimization converges (= 0 if converges)}
  \item{Wald.Ts.Longonly}{A matrix containing the Wald statistics for the testing of the longitudinal sub-model only, to examine whether the covariate is associated with the temporal presence(2nd row), or temporal non-zero proportions (3rd row), or either of the two for each covariate (1st row). Each column represents the results for one covariate. The matrix will reduce to one row for testing of the non-zero proportions if reduced model is fitted (i.e., if the zero proportion is lower than the threshold)}
  \item{pvals.Wald.Longonly}{A matrix containing the corresponding p-values of the Wald statistics from the longitudinal sub-model only}
  \item{Wald.cross.part.corr}{The Wald statistics for the testing of cross-part correlation in the longitudinal sub-model.}
  \item{pvals.Wald.Longonly}{The corresponding p-value of the Wald statistics for the testing of cross-part correlation in the longitudinal sub-model.}
}
\references{
Hu J, Wang C, Blaser M, Li H (2021). Joint modeling of zero-inflated longitudinal proportions and time-to-event data with application to a gut microbiome study. Biometrics(Accepted).
}
\author{
Jiyuan Hu
}

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
\keyword{ joint modeling }
\keyword{ microbiome }
