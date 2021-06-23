\name{dat}
\alias{dat}
\docType{data}
\title{
A demo dataset for joint modeling analysis.
}
\description{
This dataset is a data frame which contains the longitudinal microbial proportions of one taxon, the survival outcome, and (time-dependent) covariates for each subject in the long format.
}
\usage{data("dat")}
\format{
  A data frame with 467 observations on the following 9 variables. It is organized in the `long format` for both the longitudinal information and the survival information for each subject
  \describe{
    \item{\code{id}}{a factor with levels indicating the subject IDs. There are 100 subjects in this demo dataset}  
    \item{\code{Y}}{a numeric vector indicating the relative abundances (microbiobial proportions) of the investigated taxon for the subject at the observed time (obstime)}
    \item{\code{obstime}}{a numeric vector indicating the time points at which these microbiome samples were sampled}
    \item{\code{treatment}}{a numeric vector indicating the treatment group the subject is assigned to when entering the study, which is fixed througout the follow-up period for each subject. It is included as one of the covariates to adjust for}
    \item{\code{time}}{a numeric vector indicating the event time}
    \item{\code{disease}}{a numeric vector which is the indicator of the event status. It equals 1 if the event of interest occured for the corresponding subject, equals 0 if the event of interest is censored}
    \item{\code{start}}{start, stop, and event denote the risk interval and event status indicator, using the counting process formulation. The start and stop variables denote the limits of the time intervals during which the microbiome sample is collected. The start column is infact a copy of the `obstime`}
    \item{\code{stop}}{the stop variables denote the upper limits of the time intervals during which the microbiome sample is collected}
    \item{\code{event}}{a numeric vector is an indicator of whether the event occurred at the end of the corresponding time interval. It equals 1 if an event occurred at the end of the corresponding time interval}
  }
}
\details{
This demo dataset illustrates how to organize the longitudinal and survival data for the joint modeling analysis
}

\references{
Hu J, Wang C, Blaser M, Li H (2021). Joint modeling of zero-inflated longitudinal proportions and time-to-event data with application to a gut microbiome study. Biometrics(Accepted).
}
\examples{
data(dat)
## str(dat)  
}
\keyword{datasets}
