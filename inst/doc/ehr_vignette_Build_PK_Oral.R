## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----EHR_load-----------------------------------------------------------------
library(EHR)

## -----------------------------------------------------------------------------
# Data generating function for examples
mkdat <- function() {
  npat=3
  visits <- floor(runif(npat, min=2, max=6))
  id <- rep(1:npat, visits)
  dt <- as.POSIXct(paste(as.Date(sort(sample(700, sum(visits))), 
                                 origin = '2019-01-01'), '10:00:00'), tz = 'UTC') 
  + rnorm(sum(visits), 0, 1*60*60)
  dose_morn <- sample(c(2.5,5,7.5,10), sum(visits), replace = TRUE)
  conc <- round(rnorm(sum(visits), 1.5*dose_morn, 1),1)
  ld <- dt - sample(10:16, sum(visits), replace = TRUE) * 3600
  ld[rnorm(sum(visits)) < .3] <- NA
  age <- rep(sample(40:75, npat), visits)
  weight <- rep(round(rnorm(npat, 180, 20)),visits)
  hgb <- round(rep(rnorm(npat, 10, 2), visits),1)
  data.frame(id, dt, dose_morn, conc, age, weight, hgb, ld)
}

# Make example data
set.seed(30)
dat <- mkdat()
dat

## -----------------------------------------------------------------------------
dat2 <- dat[,-8]
# Build PK data without last-dose times
run_Build_PK_Oral(x = dat2,
                  idCol = "id",
                  dtCol = "dt",
                  doseCol = "dose_morn",
                  concCol = "conc",
                  ldCol = NULL,
                  first_interval_hours = 336,
                  imputeClosest = NULL)

## -----------------------------------------------------------------------------
# Build PK data with last-dose times
run_Build_PK_Oral(x = dat,
                  idCol = "id",
                  dtCol = "dt",
                  doseCol = "dose_morn",
                  concCol = "conc",
                  ldCol = "ld",
                  first_interval_hours = 336,
                  imputeClosest = NULL)


