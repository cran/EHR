## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(R.options = list(width = 100))

library(EHR)
library(lubridate)
library(pkdata)
oldopt <- options(stringsAsFactors = FALSE)

## ----niceroutput, echo = FALSE--------------------------------------------------------------------
findbreaks <- function(x, char = '[ /\\]', charlen = 75) {
  if(length(x) > 1) {
    out <- vapply(x, findbreaks, character(1), char, charlen, USE.NAMES = FALSE)
    ix <- !grepl('\n[[:space:]]*$', out)
    ix[length(ix)] <- FALSE
    out[ix] <- paste0(out[ix], '\n')
    return(paste(out, collapse = ''))
  }
  cur <- x
  nbuf <- ceiling(nchar(x) / charlen)
  if(nbuf == 1) {
    return(cur)
  }
  strings <- character(nbuf)
  i <- 1
  while(nchar(cur) > charlen) {
    loc <- c(gregexpr(char, cur)[[1]])
    b <- loc[max(which(loc < charlen))]
    strings[i] <- substr(cur, 1, b)
    cur <- substring(cur, b + 1)
    i <- i + 1
  }
  strings[i] <- cur
  paste(c(strings[1], paste0('     ', strings[-1])), collapse = '\n')
}

co <- function(expr, type = c('output', 'message')) {
  txt <- capture.output(expr, type = type)
  cat(findbreaks(txt))
  cat("\n")
}

## ----load-lib-dir, eval=FALSE---------------------------------------------------------------------
#  # load EHR package and dependencies
#  library(EHR)
#  library(pkdata)
#  library(lubridate)

## ----simp-in--------------------------------------------------------------------------------------
# define directories
td <- tempdir()
checkDir <- file.path(td, 'check1')
rawDataDir <- system.file("examples", "str_ex1", package="EHR")
dir.create(checkDir)

# pre-processed demographic data 
demo <- read.csv(file.path(rawDataDir,"Demographics_DATA_simple.csv"))
head(demo)

conc.data <- read.csv(file.path(rawDataDir,"Concentration_DATA_simple.csv"))
head(conc.data)

ivdose.data <- read.csv(file.path(rawDataDir,"IVDose_DATA_simple.csv"))
head(ivdose.data)

creat.data <- read.csv(file.path(rawDataDir,"Creatinine_DATA_simple.csv"))
head(creat.data)

## ----simp-rename----------------------------------------------------------------------------------
names(conc.data)[1:2] <- names(demo)[1:2] <- c("mod_id", "mod_id_visit")
names(creat.data)[1] <- names(ivdose.data)[1] <- "mod_id"

## ----sim-build-pk-iv, eval = FALSE----------------------------------------------------------------
#  simple_pk_dat <- run_Build_PK_IV(
#      conc=conc.data,
#      conc.columns = list(id = 'mod_id', datetime = 'date.time', druglevel = 'conc.level',
#                          idvisit = 'mod_id_visit'),
#      dose=ivdose.data,
#      dose.columns = list(id = 'mod_id', date = 'date.dose', infuseDatetime = 'infuse.time',
#                          infuseDose = 'infuse.dose', infuseTimeExact= 'infuse.time.real',
#                          bolusDatetime = 'bolus.time', bolusDose = 'bolus.dose',
#                          gap = 'maxint', weight = 'weight'),
#      demo.list = demo,
#      demo.columns = list(id = 'mod_id', idvisit = 'mod_id_visit'),
#      lab.list = list(creat.data),
#      lab.columns = list(id = 'mod_id', datetime = 'date.time'),
#      drugname='fent',
#      check.path=checkDir)

## ----sim-build-pk-iv-fake, echo = FALSE, message = FALSE------------------------------------------
co({
simple_pk_dat <- run_Build_PK_IV(
    conc=conc.data,
    conc.columns = list(id = 'mod_id', datetime = 'date.time', druglevel = 'conc.level', 
                        idvisit = 'mod_id_visit'),
    dose=ivdose.data,
    dose.columns = list(id = 'mod_id', date = 'date.dose', infuseDatetime = 'infuse.time', 
                        infuseDose = 'infuse.dose', infuseTimeExact= 'infuse.time.real', 
                        bolusDatetime = 'bolus.time', bolusDose = 'bolus.dose', 
                        gap = 'maxint', weight = 'weight'),
    demo.list = demo,
    demo.columns = list(id = 'mod_id', idvisit = 'mod_id_visit'),
    lab.list = list(creat.data),
    lab.columns = list(id = 'mod_id', datetime = 'date.time'),
    drugname='fent',
    check.path=checkDir)
}, 'message')

## ----sim-build-pk-iv-out--------------------------------------------------------------------------
head(simple_pk_dat,15)

## ----ex2-dirs-------------------------------------------------------------------------------------
dataDir <- file.path(td, 'data2')
checkDir <- file.path(td, 'check2')
rawDataDir <- system.file("examples", "str_ex2", package="EHR")
dir.create(dataDir)
dir.create(checkDir)

## ----demo-in--------------------------------------------------------------------------------------
# demographics data
demo.in <- readTransform(file.path(rawDataDir, "Demographics_DATA.csv"))
head(demo.in)

## ----samp-in1-------------------------------------------------------------------------------------
# concentration sampling times data
# read in raw data
samp.raw <- read.csv(file.path(rawDataDir, "SampleTimes_DATA.csv"))
head(samp.raw)

# transform data
samp.in0 <- dataTransformation(samp.raw,
    rename = c('Study.ID' = 'subject_id'),
    modify = list(samp = expression(as.numeric(sub('Sample ', '', Event.Name)))))
head(samp.in0)

## ----samp-in2-------------------------------------------------------------------------------------
# read in and transform data
samp.in <- readTransform(file.path(rawDataDir, "SampleTimes_DATA.csv"),
    rename = c('Study.ID' = 'subject_id'),
    modify = list(samp = expression(as.numeric(sub('Sample ', '', Event.Name)))))
head(samp.in)

## ----conc-in1-------------------------------------------------------------------------------------
# concentration sample values data
# read in raw data
conc.raw <-read.csv(file.path(rawDataDir, "SampleConcentration_DATA.csv"))
head(conc.raw)

# helper function used to make subject_id
sampId <- function(x) {
  # remove leading zeroes or trailing periods
  subid <- gsub('(^0*|\\.$)', '', x)
  # change _ to .
  gsub('_([0-9]+[_].*)$', '.\\1', subid)
}

# transform data
conc.in0 <- dataTransformation(conc.raw,
                    modify = list(
                    subid = expression(sampId(name)),
                    subject_id = expression(as.numeric(sub('[_].*', '', subid))),
                    samp = expression(sub('[^_]*[_]', '', subid)),
                    name = NULL,
                    data_file = NULL,
                    subid = NULL
                    )
                  )
head(conc.in0)

## ----conc-in2-------------------------------------------------------------------------------------
# equivalent using readTransform()
conc.in <- readTransform(file.path(rawDataDir, "SampleConcentration_DATA.csv"),
  modify = list(
    subid = expression(sampId(name)),
    subject_id = expression(as.numeric(sub('[_].*', '', subid))),
    samp = expression(sub('[^_]*[_]', '', subid)),
    name = NULL,
    data_file = NULL,
    subid = NULL
    )
  )
head(conc.in)

## ----flow-in--------------------------------------------------------------------------------------
# FLOW dosing data
flow.in <- readTransform(file.path(rawDataDir, "FLOW_DATA.csv"),
 rename = c('Subject.Id' = 'subject_id',
            'Subject.Uniq.Id' = 'subject_uid'),
 modify=list(
  date.time = expression(pkdata::parse_dates(EHR:::fixDates(Perform.Date))),
  unit = expression(sub('.*[ ]', '', Final.Rate..NFR.units.)),
  rate = expression(as.numeric(sub('([0-9.]+).*', '\\1', Final.Rate..NFR.units.)))
  )
 ) 
head(flow.in)

## ----flow-in0-------------------------------------------------------------------------------------
# FLOW dosing data
flow.in <- readTransform(file.path(rawDataDir, "FLOW_DATA.csv"),
                         rename = c('Subject.Id' = 'subject_id',
                                    'Subject.Uniq.Id' = 'subject_uid')) 
# pre-process the flow data 
# date.time variable should be in an appropriate form
flow.in[,'date.time'] <- pkdata::parse_dates(EHR:::fixDates(flow.in[,'Perform.Date']))
# unit and rate are required: separate unit and rate from 'Final.Rate..NFR.units.' if needed
flow.in[,'unit'] <- sub('.*[ ]', '', flow.in[,'Final.Rate..NFR.units.'])
flow.in[,'rate'] <- as.numeric(sub('([0-9.]+).*', '\\1', flow.in[,'Final.Rate..NFR.units.']))
head(flow.in)

## ----mar-in---------------------------------------------------------------------------------------
# MAR dosing data
mar.in0 <- read.csv(file.path(rawDataDir, "MAR_DATA.csv"), check.names = FALSE)
mar.in <- dataTransformation(mar.in0, rename = c('Uniq.Id' = 'subject_uid'))
head(mar.in)

## ----labs-in--------------------------------------------------------------------------------------
# Serum creatinine lab data
creat.in <- readTransform(file.path(rawDataDir, "Creatinine_DATA.csv"),
    rename = c('Subject.uniq' = 'subject_uid'))
head(creat.in)

# Albumin lab data
alb.in <- readTransform(file.path(rawDataDir, "Albumin_DATA.csv"),
    rename = c('Subject.uniq' = 'subject_uid'))
head(alb.in)

## ----merge-ids------------------------------------------------------------------------------------
# merge all ID datasets
data <-  list(demo.in,
              samp.in,
              conc.in,
              flow.in,
              mar.in,
              creat.in,
              alb.in)

idcols <-  list(c('subject_id', 'subject_uid'), # id vars in demo.in
                'subject_id', # id var in samp.in
                'subject_id', # id var in conc.in
                c('subject_id', 'subject_uid'), # id vars in flow.in
                'subject_uid', # id var in mar.in
                'subject_uid', # id var in creat.in
                'subject_uid') # id var in creat.in

mod.id <- idCrosswalk(data, idcols, visit.id="subject_id", uniq.id="subject_uid")
saveRDS(mod.id, file=file.path(dataDir,"Fentanyl_module_id.rds"))

mod.id

## ---- eval = FALSE--------------------------------------------------------------------------------
#  pullFakeId(dat, xwalk, firstCols = NULL, orderBy = NULL)

## ----mod-id-data----------------------------------------------------------------------------------
## demographics data
demo.cln <- pullFakeId(demo.in, mod.id,
    firstCols = c('mod_id', 'mod_visit', 'mod_id_visit'),
    uniq.id = 'subject_uid')
head(demo.cln)
saveRDS(demo.cln, file=file.path(dataDir,"Fentanyl_demo_mod_id.rds"))

## drug level data
# sampling times
samp.cln <- pullFakeId(samp.in, mod.id,
    firstCols = c('mod_id', 'mod_visit', 'mod_id_visit', 'samp'), 
    orderBy = c('mod_id_visit','samp'),
    uniq.id = 'subject_uid')
head(samp.cln)
saveRDS(samp.cln, file=file.path(dataDir,"Fentanyl_samp_mod_id.rds"))

# sampling concentrations
conc.cln <- pullFakeId(conc.in, mod.id,
    firstCols = c('record_id', 'mod_id', 'mod_visit', 'mod_id_visit', 'samp'),
    orderBy = 'record_id',
    uniq.id = 'subject_uid')
head(conc.cln)
saveRDS(conc.cln, file=file.path(dataDir,"Fentanyl_conc_mod_id.rds"))

## dosing data
# flow
flow.cln <- pullFakeId(flow.in, mod.id,
    firstCols = c('mod_id', 'mod_visit', 'mod_id_visit'),
    uniq.id = 'subject_uid')
head(flow.cln)
saveRDS(flow.cln, file=file.path(dataDir,"Fentanyl_flow_mod_id.rds"))

# mar
mar.cln <- pullFakeId(mar.in, mod.id, firstCols = 'mod_id', uniq.id = 'subject_uid')
head(mar.cln)
saveRDS(mar.cln, file=file.path(dataDir,"Fentanyl_mar_mod_id.rds"))

## laboratory data
creat.cln <- pullFakeId(creat.in, mod.id, 'mod_id',uniq.id = 'subject_uid')
head(creat.cln)

alb.cln <- pullFakeId(alb.in, mod.id, 'mod_id', uniq.id = 'subject_uid')
head(alb.cln)

saveRDS(creat.cln, file=file.path(dataDir,"Fentanyl_creat_mod_id.rds"))
saveRDS(alb.cln, file=file.path(dataDir,"Fentanyl_alb_mod_id.rds"))

## ----mod-setup------------------------------------------------------------------------------------
# set crosswalk option 
xwalk <- readRDS(file.path(dataDir, "Fentanyl_module_id.rds"))
oldxwalk <- options(pkxwalk = 'xwalk')

# define parameters
drugname <- 'fent'
LLOQ <- 0.05

## ----Pro-Demographic, eval=TRUE-------------------------------------------------------------------
# helper function
exclude_val <- function(x, val=1) { !is.na(x) & x == val }

demo.out <- run_Demo(demo.path = file.path(dataDir, "Fentanyl_demo_mod_id.rds"),
    demo.columns = list(id = 'mod_id_visit'),
    toexclude = expression(exclude_val(in_hospital_mortality) | exclude_val(add_ecmo)),
    demo.mod.list = list(length_of_icu_stay = 
                        expression(daysDiff(surgery_date, date_icu_dc))))

head(demo.out$demo)
demo.out$exclude

## ----Pro-Med-Str1, eval=FALSE---------------------------------------------------------------------
#  ivdose.out <- run_MedStrI(
#      mar.path=file.path(dataDir,"Fentanyl_mar_mod_id.rds"),
#      mar.columns = list(id = 'mod_id', datetime = c('Date','Time'), dose = 'med:dosage',
#                         drug = 'med:mDrug', given = 'med:given'),
#      medGivenReq = TRUE,
#      flow.path=file.path(dataDir,"Fentanyl_flow_mod_id.rds"),
#      flow.columns = list(id = 'mod_id', datetime = 'date.time', finalunits = 'Final.Units',
#                          unit = 'unit', rate = 'rate', weight = 'Final.Wt..kg.'),
#      medchk.path = file.path(rawDataDir, sprintf('medChecked-%s.csv', drugname)),
#      demo.list = NULL,
#      demo.columns = list(),
#      missing.wgt.path = NULL,
#      wgt.columns = list(),
#      check.path = checkDir,
#      failflow_fn = 'FailFlow',
#      failunit_fn = 'Unit',
#      failnowgt_fn = 'NoWgt',
#      infusion.unit = 'mcg/kg/hr',
#      bolus.unit = 'mcg',
#      bol.rate.thresh = Inf,
#      rateunit = 'mcg/hr',
#      ratewgtunit = 'mcg/kg/hr',
#      weightunit = 'kg',
#      drugname = drugname)

## ----Pro-Med-Str1-fake, echo=FALSE, message = FALSE-----------------------------------------------
co({
ivdose.out <- run_MedStrI(
    mar.path=file.path(dataDir,"Fentanyl_mar_mod_id.rds"),
    mar.columns = list(id = 'mod_id', datetime = c('Date','Time'), dose = 'med:dosage', drug = 'med:mDrug', given = 'med:given'),
    medGivenReq = TRUE,
    flow.path=file.path(dataDir,"Fentanyl_flow_mod_id.rds"),
    flow.columns = list(id = 'mod_id', datetime = 'date.time', finalunits = 'Final.Units', unit = 'unit', rate = 'rate', weight = 'Final.Wt..kg.'),
    medchk.path = file.path(rawDataDir, sprintf('medChecked-%s.csv', drugname)),
    demo.list = NULL,
    demo.columns = list(),
    missing.wgt.path = NULL,
    wgt.columns = list(),
    check.path = checkDir,
    failflow_fn = 'FailFlow',
    failunit_fn = 'Unit',
    failnowgt_fn = 'NoWgt',
    infusion.unit = 'mcg/kg/hr',
    bolus.unit = 'mcg',
    bol.rate.thresh = Inf,
    rateunit = 'mcg/hr',
    ratewgtunit = 'mcg/kg/hr',
    weightunit = 'kg',
    drugname = drugname)
}, 'message')

## ----Pro-Med-Str1-out-----------------------------------------------------------------------------
head(ivdose.out)

## ----eRX-dat--------------------------------------------------------------------------------------
(eRX <- read.csv(file.path(rawDataDir,"e-rx_DATA.csv"),stringsAsFactors = FALSE))

## ----Pro-Med-Str2---------------------------------------------------------------------------------
eRX.out <- run_MedStrII(file.path(rawDataDir,"e-rx_DATA.csv"),
    dat.columns = list(id = 'GRID', dose = 'RX_DOSE', freq = 'FREQUENCY', date = 'ENTRY_DATE', 
                       str = 'STRENGTH_AMOUNT', desc = 'DESCRIPTION')
)

eRX.out

## ----Pro-Drug-Level, eval=FALSE-------------------------------------------------------------------
#  conc.out <- run_DrugLevel(conc.path=file.path(dataDir,"Fentanyl_conc_mod_id.rds"),
#      conc.columns = list(
#        id = 'mod_id', conc = 'conc.level',
#        idvisit = 'mod_id_visit', samplinkid = 'mod_id_event'
#      ),
#      conc.select=c('mod_id','mod_id_visit','samp','fentanyl_calc_conc'),
#      conc.rename=c(fentanyl_calc_conc = 'conc.level', samp= 'event'),
#      conc.mod.list=list(mod_id_event = expression(paste(mod_id_visit, event, sep = '_'))),
#      samp.path=file.path(dataDir,"Fentanyl_samp_mod_id.rds"),
#      samp.columns = list(
#        conclinkid = 'mod_id_event', datetime = 'Sample.Collection.Date.and.Time'
#      ),
#      samp.mod.list=list(mod_id_event = expression(paste(mod_id_visit, samp, sep = '_'))),
#      check.path=checkDir,
#      failmiss_fn = 'MissingConcDate-',
#      multsets_fn = 'multipleSetsConc-',
#      faildup_fn = 'DuplicateConc-',
#      drugname = drugname,
#      LLOQ = LLOQ,
#      demo.list = demo.out,
#      demo.columns = list(id = 'mod_id', idvisit = 'mod_id_visit'))

## ----Pro-Drug-Level-fake, echo=FALSE--------------------------------------------------------------
co({
conc.out <- run_DrugLevel(conc.path=file.path(dataDir,"Fentanyl_conc_mod_id.rds"),
    conc.columns = list(id = 'mod_id', conc = 'conc.level', idvisit = 'mod_id_visit', samplinkid = 'mod_id_event'),
    conc.select=c('mod_id','mod_id_visit','samp','fentanyl_calc_conc'),
    conc.rename=c(fentanyl_calc_conc = 'conc.level', samp= 'event'),
    conc.mod.list=list(mod_id_event = expression(paste(mod_id_visit, event, sep = '_'))),
    samp.path=file.path(dataDir,"Fentanyl_samp_mod_id.rds"),
    samp.columns = list(conclinkid = 'mod_id_event', datetime = 'Sample.Collection.Date.and.Time'),
    samp.mod.list=list(mod_id_event = expression(paste(mod_id_visit, samp, sep = '_'))),
    check.path=checkDir,
    failmiss_fn = 'MissingConcDate-',
    multsets_fn = 'multipleSetsConc-',
    faildup_fn = 'DuplicateConc-',
    drugname = drugname,
    LLOQ = LLOQ,
    demo.list = demo.out,
    demo.columns = list(id = 'mod_id', idvisit = 'mod_id_visit'))
})

## ----Pro-Drug-Level-out---------------------------------------------------------------------------
head(conc.out)

## ----faildate-------------------------------------------------------------------------------------
( fail.miss.conc.date <- read.csv(file.path(checkDir,"failMissingConcDate-fent.csv")) )

## ----fixdate--------------------------------------------------------------------------------------
fail.miss.conc.date[,"datetime"] <- c("9/30/2016 09:32","10/1/2016 19:20","10/2/2016 02:04")
fail.miss.conc.date
 
write.csv(fail.miss.conc.date, file.path(checkDir,"fixMissingConcDate-fent.csv"))

## ----Pro-Drug-Level-rerun, eval=FALSE-------------------------------------------------------------
#  conc.out <- run_DrugLevel(conc.path=file.path(dataDir,"Fentanyl_conc_mod_id.rds"),
#      conc.columns = list(
#        id = 'mod_id', conc = 'conc.level',
#        idvisit = 'mod_id_visit', samplinkid = 'mod_id_event'
#      ),
#      conc.select=c('mod_id','mod_id_visit','samp','fentanyl_calc_conc'),
#      conc.rename=c(fentanyl_calc_conc = 'conc.level', samp = 'event'),
#      conc.mod.list=list(mod_id_event = expression(paste(mod_id_visit, event, sep = '_'))),
#      samp.path=file.path(dataDir,"Fentanyl_samp_mod_id.rds"),
#      samp.columns = list(
#        conclinkid = 'mod_id_event', datetime = 'Sample.Collection.Date.and.Time'
#      ),
#      samp.mod.list=list(mod_id_event = expression(paste(mod_id_visit, samp, sep = '_'))),
#      check.path=checkDir,
#      failmiss_fn = 'MissingConcDate-',
#      multsets_fn = 'multipleSetsConc-',
#      faildup_fn = 'DuplicateConc-',
#      drugname = drugname,
#      LLOQ = LLOQ,
#      demo.list = demo.out,
#      demo.columns = list(id = 'mod_id', idvisit = 'mod_id_visit'))

## ----Pro-Drug-Level-rerun-fake, echo=FALSE--------------------------------------------------------
co({
conc.out <- run_DrugLevel(conc.path=file.path(dataDir,"Fentanyl_conc_mod_id.rds"),
    conc.columns = list(id = 'mod_id', conc = 'conc.level', idvisit = 'mod_id_visit', samplinkid = 'mod_id_event'),
    conc.select=c('mod_id','mod_id_visit','samp','fentanyl_calc_conc'),
    conc.rename=c(fentanyl_calc_conc = 'conc.level', samp = 'event'),
    conc.mod.list=list(mod_id_event = expression(paste(mod_id_visit, event, sep = '_'))),
    samp.path=file.path(dataDir,"Fentanyl_samp_mod_id.rds"),
    samp.columns = list(conclinkid = 'mod_id_event', datetime = 'Sample.Collection.Date.and.Time'),
    samp.mod.list=list(mod_id_event = expression(paste(mod_id_visit, samp, sep = '_'))),
    check.path=checkDir,
    failmiss_fn = 'MissingConcDate-',
    multsets_fn = 'multipleSetsConc-',
    faildup_fn = 'DuplicateConc-', 
    drugname = drugname,
    LLOQ = LLOQ,
    demo.list = demo.out,
    demo.columns = list(id = 'mod_id', idvisit = 'mod_id_visit'))
})

## ----remove-fix, include=FALSE--------------------------------------------------------------------
# remove fix file, so running vignette produces warning with first run of run_DrugLevel()
fx <- file.path(checkDir,"fixMissingConcDate-fent.csv")
if (file.exists(fx)) file.remove(fx)

# remove multiplesetsconc file
ms <- file.path(checkDir,paste0("multipleSetsConc-", drugname, Sys.Date(),".csv"))
if (file.exists(ms)) file.remove(ms)

## ----Pro-Laboratory, eval=TRUE--------------------------------------------------------------------
creat.out <- run_Labs(lab.path=file.path(dataDir,"Fentanyl_creat_mod_id.rds"),
    lab.select = c('mod_id','date.time','creat'),
    lab.mod.list = list(date.time = expression(parse_dates(fixDates(paste(date, time))))))

alb.out <- run_Labs(lab.path=file.path(dataDir,"Fentanyl_alb_mod_id.rds"),
    lab.select = c('mod_id','date.time','alb'),
    lab.mod.list = list(date.time = expression(parse_dates(fixDates(paste(date, time))))))

lab.out <- list(creat.out, alb.out)

str(lab.out)

## ----Build-PK-IV, eval=FALSE----------------------------------------------------------------------
#  pk_dat <- run_Build_PK_IV(
#      conc=conc.out,
#      conc.columns = list(id = 'mod_id', datetime = 'date.time', druglevel = 'conc.level',
#                          idvisit = 'mod_id_visit'),
#      dose=ivdose.out,
#      dose.columns = list(id = 'mod_id', date = 'date.dose', infuseDatetime = 'infuse.time',
#                          infuseDose = 'infuse.dose', infuseTimeExact= 'infuse.time.real',
#                          bolusDatetime = 'bolus.time', bolusDose = 'bolus.dose',
#                          gap = 'maxint', weight = 'weight'),
#      demo.list = demo.out,
#      demo.columns = list(id = 'mod_id', idvisit = 'mod_id_visit'),
#      lab.list = lab.out,
#      lab.columns = list(id = 'mod_id', datetime = 'date.time'),
#      pk.vars=c('date'),
#      drugname=drugname,
#      check.path=checkDir,
#      missdemo_fn='-missing-demo',
#      faildupbol_fn='DuplicateBolus-',
#      date.format="%m/%d/%y %H:%M:%S",
#      date.tz="America/Chicago")

## ----Build-PK-IV-fake, echo=FALSE, message=FALSE--------------------------------------------------
co({
pk_dat <- run_Build_PK_IV(
    conc=conc.out,
    conc.columns = list(id = 'mod_id', datetime = 'date.time', druglevel = 'conc.level', 
                        idvisit = 'mod_id_visit'),
    dose=ivdose.out,
    dose.columns = list(id = 'mod_id', date = 'date.dose', infuseDatetime = 'infuse.time', 
                        infuseDose = 'infuse.dose', infuseTimeExact= 'infuse.time.real',  bolusDatetime = 'bolus.time', bolusDose = 'bolus.dose', gap = 'maxint', weight = 'weight'),
    demo.list = demo.out,
    demo.columns = list(id = 'mod_id', idvisit = 'mod_id_visit'),
    lab.list = lab.out,
    lab.columns = list(id = 'mod_id', datetime = 'date.time'),
    pk.vars=c('date'),
    drugname=drugname,
    check.path=checkDir,
    missdemo_fn='-missing-demo',
    faildupbol_fn='DuplicateBolus-',
    date.format="%m/%d/%y %H:%M:%S",
    date.tz="America/Chicago")
}, 'message')

## ----Build-PK-IV-out------------------------------------------------------------------------------
# convert id back to original IDs
pk_dat <- pullRealId(pk_dat, remove.mod.id=TRUE)

head(pk_dat)

## ----cleanup, echo = FALSE------------------------------------------------------------------------
# restore old options
options(oldopt)
options(oldxwalk)
# normally, you would not delete these files
# CRAN policy states that a package should do proper cleanup
to_delete <- c(file.path(td, 'check1'), file.path(td, 'check2'), file.path(td, 'data2'))
unlink(to_delete, recursive = TRUE)

