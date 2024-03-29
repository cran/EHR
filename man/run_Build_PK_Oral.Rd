% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkDosingInterval.R
\name{run_Build_PK_Oral}
\alias{run_Build_PK_Oral}
\title{Build-PK-Oral Module}
\usage{
run_Build_PK_Oral(
  x,
  idCol = "id",
  dtCol = "dt",
  doseCol = "dose",
  concCol = "conc",
  ldCol = NULL,
  first_interval_hours = 336,
  imputeClosest = NULL
)
}
\arguments{
\item{x}{a data.frame or file saved as either CSV, RData, or RDS}

\item{idCol}{data.frame id column name}

\item{dtCol}{data.frame date column name}

\item{doseCol}{dose column name}

\item{concCol}{concentration column name}

\item{ldCol}{last-dose time column name}

\item{first_interval_hours}{number of hours before the first concentration to start time=0; the default is 336 hours = 14 days}

\item{imputeClosest}{columns to impute missing data with next observation propagated backward; this is in addition to
all covariates receving imputation using last observation carried forward}
}
\value{
data.frame
}
\description{
This module builds PK data for orally administered medications.
}
\details{
See EHR Vignette for Build-PK-Oral.
}
\examples{
## Data Generating Function
mkdat <- function() {
  npat <- 3
  visits <- floor(runif(npat, min=2, max=6))
  id <- rep(1:npat, visits)
  dt_samp <- as.Date(sort(sample(700, sum(visits))), origin = '2019-01-01')
  tm_samp <- as.POSIXct(paste(dt_samp, '10:00:00'), tz = 'UTC')
  dt <- tm_samp + rnorm(sum(visits), 0, 1*60*60)
  dose_morn <- sample(c(2.5,5,7.5,10), sum(visits), replace = TRUE)
  conc <- round(rnorm(sum(visits), 1.5*dose_morn, 1),1)
  ld <- dt - sample(10:16, sum(visits), replace = TRUE) * 3600
  ld[rnorm(sum(visits)) < .3] <- NA
  age <- rep(sample(40:75, npat), visits)
  gender <- rep(sample(0:1, npat, replace=TRUE), visits)
  weight <- rep(round(rnorm(npat, 180, 20)),visits)
  hgb <- rep(rnorm(npat, 10, 2), visits)
  data.frame(id, dt, dose_morn, conc, ld, age, gender, weight, hgb)
}

# Make raw data
set.seed(30)
dat <- mkdat()

#Process data without last-dose times
run_Build_PK_Oral(x = dat,
                  idCol = "id",
                  dtCol = "dt",
                  doseCol = "dose_morn",
                  concCol = "conc",
                  ldCol = NULL,
                  first_interval_hours = 336,
                  imputeClosest = NULL)

#Process data with last-dose times
run_Build_PK_Oral(x = dat, doseCol = "dose_morn", ldCol = "ld")

}
