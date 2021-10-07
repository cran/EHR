## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(EHR)

## ----niceroutput, echo = FALSE------------------------------------------------
breaktable <- function(df, brks) {
  if(class(df)[1] != 'data.frame') {
    df <- as.data.frame(df)
  }
  if(max(brks) != ncol(df)) {
    brks <- sort(c(brks, ncol(df)))
  }
  lb <- length(brks)
  res <- vector('list', lb * 2 - 1)
  pos <- 1
  for(i in seq(lb)) {
    curdf <- df[, seq(pos, brks[i]), drop = FALSE]
    curout <- capture.output(print(curdf))
    if(i > 1) curout <- paste0('     ', curout)
    res[[(i - 1) * 2 + 1]] <- curout
    if(i < lb) res[[i * 2]] <- ''
    pos <- brks[i] + 1
  }
  paste(do.call(c, res), collapse = '\n')
}

## -----------------------------------------------------------------------------
tac_mxr_fn <- system.file("examples", "tac_mxr_out.csv", package = "EHR")
tac_mxr <- read.csv(tac_mxr_fn, na = '')
tac_mxr[c(135:139,163:167,283:289,343:346),]

## -----------------------------------------------------------------------------
tac_mxr_parsed <- parseMedExtractR(tac_mxr_fn)

## ----echo=FALSE---------------------------------------------------------------
cat(breaktable(tac_mxr_parsed[c(32,33,42,43,87,88,102),], c(3,7)))

## -----------------------------------------------------------------------------
tac_mxr_part1_out <- buildDose(tac_mxr_parsed)

## ----echo=FALSE---------------------------------------------------------------
tac_mxr_part1_out[c(51,52,55,56,104,105,106),]

## -----------------------------------------------------------------------------
tac_gs_part1 <- read.csv(system.file("examples", "tac_gs_part1.csv", package = "EHR"),
                         stringsAsFactors = FALSE, na = '')

## ----echo=FALSE---------------------------------------------------------------
tac_gs_part1[c(51:54,104,105,107),]

## -----------------------------------------------------------------------------
precall <- function(dat, gs) {
  tp1 <- sum(dat %in% gs)
  fp1 <- sum(!(dat %in% gs))
  fn1 <- sum(!(gs %in% dat))
  r1 <- c(tp1, tp1 + fn1)
  p1 <- c(tp1, tp1 + fp1)
  r <- rbind(r1,p1)
  dimnames(r) <- list(c('recall','prec'), c('num','den'))
  cbind(r, prop = round(r[,1] / r[,2], 2))
}

colsToCompare <- c('filename','drugname','strength','dose','route','freq',
  'dosestr','dosechange','drugname_start')
tac_mxr_part1_out <- tac_mxr_part1_out[,colsToCompare]
tac_gs_part1 <- tac_gs_part1[,colsToCompare]

tacxrrow <- do.call(paste, c(tac_mxr_part1_out, sep = '|'))
gs.tacxrrow <- do.call(paste, c(tac_gs_part1, sep = '|'))

precall(tacxrrow, gs.tacxrrow)

## -----------------------------------------------------------------------------
bmd <- function(x) {
  fns <- strsplit(x, '_')
  pid <- sapply(fns, `[`, 1)
  date <- as.Date(sapply(fns, `[`, 2), format = '%Y-%m-%d')
  note <- sapply(fns, `[`, 3)
  data.frame(filename = x, pid, date, note, stringsAsFactors = FALSE)
}
tac_metadata <- bmd(tac_mxr_part1_out[['filename']])

## ----echo=FALSE---------------------------------------------------------------
tac_metadata[c(51,55,104,105),]

## ---- eval = FALSE------------------------------------------------------------
#  tac_part2 <- collapseDose(tac_mxr_part1_out, tac_metadata, naFreq='most')

## ---- echo = FALSE, message = FALSE-------------------------------------------
suppressWarnings(tac_part2 <- collapseDose(tac_mxr_part1_out, tac_metadata, naFreq='most'))

## ----echo=FALSE---------------------------------------------------------------
cat(breaktable(tac_part2$note[c(40,42,68,69),], c(7,12)))

## ----echo=FALSE---------------------------------------------------------------
cat(breaktable(tac_part2$date[c(29,42),], c(7,12)))

## -----------------------------------------------------------------------------
tac_gs_part2_note <- read.csv(
  system.file("examples", "tac_gs_part2_note.csv", package = "EHR"),
  stringsAsFactors = FALSE, na = ''
)

## ----echo=FALSE---------------------------------------------------------------
tac_gs_part2_note[c(40,41,68,70),]

## -----------------------------------------------------------------------------
tac_gs_part2_date <- read.csv(
  system.file("examples", "tac_gs_part2_date.csv", package = "EHR"),
  stringsAsFactors = FALSE, na = ''
)

## ----echo=FALSE---------------------------------------------------------------
tac_gs_part2_date[c(29,42),]

## ---- eval = FALSE------------------------------------------------------------
#  precall <- function(dat, gs) {
#    tp1 <- sum(dat %in% gs)
#    fp1 <- sum(!(dat %in% gs))
#    fn1 <- sum(!(gs %in% dat))
#    r1 <- c(tp1, tp1 + fn1)
#    p1 <- c(tp1, tp1 + fp1)
#    r <- rbind(r1,p1)
#    dimnames(r) <- list(c('recall','prec'), c('num','den'))
#    cbind(r, prop = round(r[,1] / r[,2], 2))
#  }
#  
#  metaData <- bmd(unique(tac_gs_part1$filename))
#  tacxr <- collapseDose(tac_gs_part1, metaData, 'bid')
#  tacxr.note <- tacxr[['note']]
#  tacxr.date <- tacxr[['date']]
#  
#  tacxr.note$pid <- sub("_.*","",tacxr.note$filename)
#  tacxr.date$pid <- sub("_.*","",tacxr.date$filename)
#  tac_gs_part2_note$pid <- sub("_.*","",tac_gs_part2_note$filename)
#  tac_gs_part2_date$pid <- sub("_.*","",tac_gs_part2_date$filename)
#  
#  tacxrrow.note.intake <- do.call(paste, c(tacxr.note[,c('pid','dose.intake',
#                                                         'dosechange')],sep = '|'))
#  tacxrrow.note.daily <- do.call(paste, c(tacxr.note[,c('pid','intaketime','dose.daily',
#                                                        'dosechange')], sep = '|'))
#  tacxrrow.date.intake <- do.call(paste, c(tacxr.date[,c('pid','dose.intake',
#                                                         'dosechange')], sep = '|'))
#  tacxrrow.date.daily <- do.call(paste, c(tacxr.date[,c('pid','intaketime','dose.daily',
#                                                        'dosechange')], sep = '|'))
#  
#  gs.tacxrrow.note.intake <- do.call(paste, c(tac_gs_part2_note[,c('pid','doseintake',
#                                                                   'dosechange')], sep = '|'))
#  gs.tacxrrow.note.daily <- do.call(paste, c(tac_gs_part2_note[,c('pid','intaketime','daily',
#                                                                  'dosechange')], sep = '|'))
#  gs.tacxrrow.date.intake <- do.call(paste, c(tac_gs_part2_date[,c('pid','doseintake',
#                                                                   'dosechange')], sep = '|'))
#  gs.tacxrrow.date.daily <- do.call(paste, c(tac_gs_part2_date[,c('pid','intaketime','daily',
#                                                                  'dosechange')], sep = '|'))
#  
#  precall(tacxrrow.note.intake, gs.tacxrrow.note.intake)
#  precall(tacxrrow.note.daily, gs.tacxrrow.note.daily)
#  precall(tacxrrow.date.intake, gs.tacxrrow.date.intake)
#  precall(tacxrrow.date.daily, gs.tacxrrow.date.daily)

## ---- echo = FALSE, warning = FALSE-------------------------------------------
precall <- function(dat, gs) {
  tp1 <- sum(dat %in% gs)
  fp1 <- sum(!(dat %in% gs))
  fn1 <- sum(!(gs %in% dat))
  r1 <- c(tp1, tp1 + fn1)
  p1 <- c(tp1, tp1 + fp1)
  r <- rbind(r1,p1)
  dimnames(r) <- list(c('recall','prec'), c('num','den'))
  cbind(r, prop = round(r[,1] / r[,2], 2))
}

metaData <- bmd(unique(tac_gs_part1$filename))
suppressWarnings(tacxr <- collapseDose(tac_gs_part1, metaData, 'bid'))
tacxr.note <- tacxr[['note']]
tacxr.date <- tacxr[['date']]

tacxr.note$pid <- sub("_.*","",tacxr.note$filename)
tacxr.date$pid <- sub("_.*","",tacxr.date$filename)
tac_gs_part2_note$pid <- sub("_.*","",tac_gs_part2_note$filename)
tac_gs_part2_date$pid <- sub("_.*","",tac_gs_part2_date$filename)

tacxrrow.note.intake <- do.call(paste, c(tacxr.note[,c('pid','dose.intake',
                                                       'dosechange')],sep = '|'))
tacxrrow.note.daily <- do.call(paste, c(tacxr.note[,c('pid','intaketime','dose.daily',
                                                      'dosechange')], sep = '|'))
tacxrrow.date.intake <- do.call(paste, c(tacxr.date[,c('pid','dose.intake',
                                                       'dosechange')], sep = '|'))
tacxrrow.date.daily <- do.call(paste, c(tacxr.date[,c('pid','intaketime','dose.daily',
                                                      'dosechange')], sep = '|'))

gs.tacxrrow.note.intake <- do.call(paste, c(tac_gs_part2_note[,c('pid','doseintake',
                                                                 'dosechange')], sep = '|'))
gs.tacxrrow.note.daily <- do.call(paste, c(tac_gs_part2_note[,c('pid','intaketime','daily',
                                                                'dosechange')], sep = '|'))
gs.tacxrrow.date.intake <- do.call(paste, c(tac_gs_part2_date[,c('pid','doseintake',
                                                                 'dosechange')], sep = '|'))
gs.tacxrrow.date.daily <- do.call(paste, c(tac_gs_part2_date[,c('pid','intaketime','daily',
                                                                'dosechange')], sep = '|'))

precall(tacxrrow.note.intake, gs.tacxrrow.note.intake)
precall(tacxrrow.note.daily, gs.tacxrrow.note.daily)
precall(tacxrrow.date.intake, gs.tacxrrow.date.intake)
precall(tacxrrow.date.daily, gs.tacxrrow.date.daily)

