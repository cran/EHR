## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)

## ----niceroutput, echo = FALSE------------------------------------------------
findbreaks <- function(x, char = '\\|', charlen = 75) {
  if(length(x) > 1) {
    out <- vapply(x, findbreaks, character(1), char, charlen, USE.NAMES = FALSE)
    return(paste(out, collapse = '\n'))
  }
  cur <- x
  nbuf <- ceiling(nchar(x) / charlen)
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

## ----EHR_load-----------------------------------------------------------------
library(EHR)

## ---- eval = FALSE------------------------------------------------------------
#  extractMed(note_fn, drugnames, drgunit, windowlength, max_edit_dist, ...)

## ----mxr----------------------------------------------------------------------
# tacrolimus note file names
tac_fn <- list(
  system.file("examples", "tacpid1_2008-06-26_note1_1.txt", package = "EHR"),
  system.file("examples", "tacpid1_2008-06-26_note2_1.txt", package = "EHR"),
  system.file("examples", "tacpid1_2008-12-16_note3_1.txt", package = "EHR")
)

# execute module
tac_mxr <- extractMed(tac_fn,
                       drugnames = c("tacrolimus", "prograf", "tac", "tacro", "fk", "fk506"),
                       drgunit = "mg",
                       windowlength = 60,
                       max_edit_dist = 2,
                       lastdose=TRUE)


# lamotrigine note file name
lam_fn <- c(
  system.file("examples", "lampid1_2016-02-05_note4_1.txt", package = "EHR"),
  system.file("examples", "lampid1_2016-02-05_note5_1.txt", package = "EHR"),
  system.file("examples", "lampid2_2008-07-20_note6_1.txt", package = "EHR"),
  system.file("examples", "lampid2_2012-04-15_note7_1.txt", package = "EHR")
)

# execute module
lam_mxr <- extractMed(lam_fn,
                       drugnames = c("lamotrigine", "lamotrigine XR", 
                                     "lamictal", "lamictal XR", 
                                     "LTG", "LTG XR"),
                       drgunit = "mg",
                       windowlength = 130,
                       max_edit_dist = 1,
                       strength_sep="-")

## ---- echo = FALSE------------------------------------------------------------
# Print output
message("tacrolimus medExtractR output:\n")
head(tac_mxr)

message("lamotrigine medExtractR output:\n")
head(lam_mxr)

## ---- eval = FALSE------------------------------------------------------------
#  # save as csv files
#  write.csv(tac_mxr, file='tac_mxr.csv', row.names=FALSE)
#  write.csv(lam_mxr, file='lam_mxr.csv', row.names=FALSE)

## ---- echo = FALSE------------------------------------------------------------
ann <- read.delim(system.file("mxr_tune", "ann_example.ann", package = "EHR"), 
                            header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                            col.names = c("id", "entity", "annotation"))
head(ann)

## -----------------------------------------------------------------------------
# Read in the annotations - might be specific to annotation method/software
ann_filenames <- list(system.file("mxr_tune", "tune_note1.ann", package = "EHR"),
                      system.file("mxr_tune", "tune_note2.ann", package = "EHR"),
                      system.file("mxr_tune", "tune_note3.ann", package = "EHR"))

tune_ann <- do.call(rbind, lapply(ann_filenames, function(fn){
  annotations <- read.delim(fn, 
                            header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                            col.names = c("id", "entity", "annotation"))
  
  # Label with file name
  annotations$filename <- sub(".ann", ".txt", sub(".+/", "", fn), fixed=TRUE)
  
  # Separate entity information into entity label and start:stop position
  # Format is "entity start stop"
  ent_info <- strsplit(as.character(annotations$entity), split="\\s")
  annotations$entity <- unlist(lapply(ent_info, '[[', 1))
  annotations$pos <- paste(lapply(ent_info, '[[', 2), 
                           lapply(ent_info, '[[', 3), sep=":")
  
  annotations <- annotations[,c("filename", "entity", "annotation", "pos")]
  
  return(annotations)
}))

head(tune_ann)

## ----run_mxr------------------------------------------------------------------
wind_len <- seq(30, 120, 30)
max_edit <- seq(0, 2, 1)
tune_pick <- expand.grid("window_length" = wind_len, 
                         "max_edit_distance" = max_edit)

# Run the Extract-Med module on the tuning notes
note_filenames <- list(system.file("mxr_tune", "tune_note1.txt", package = "EHR"),
                       system.file("mxr_tune", "tune_note2.txt", package = "EHR"),
                       system.file("mxr_tune", "tune_note3.txt", package = "EHR"))

# List to store output for each parameter combination
mxr_tune <- vector(mode="list", length=nrow(tune_pick))

for(i in 1:nrow(tune_pick)){
  mxr_tune[[i]] <- extractMed(note_filenames,
               drugnames = c("tacrolimus", "prograf", "tac", "tacro", "fk", "fk506"),
               drgunit = "mg",
               windowlength = tune_pick$window_length[i],
               max_edit_dist = tune_pick$max_edit_distance[i],
               progress = FALSE)
}


## -----------------------------------------------------------------------------
# Functions to compute true positive, false positive, and false negatives
# number of true positives - how many annotations were correctly identified by extractMed
Tpos <- function(df){
  sum(df$annotation == df$expr, na.rm=TRUE)
}
# number of false positive (identified by extractMed but not annotated)
Fpos <- function(df){
  sum(is.na(df$annotation))
}
# number of false negatives (annotated but not identified by extractMed)
Fneg <- function(df){
  # keep only rows with annotation
  df_ann <- subset(df, !is.na(annotation))
  sum(is.na(df$expr))
}

prf <- function(df){
  tp <- Tpos(df)
  fp <- Fpos(df)
  fn <- Fneg(df)
  
  precision <- tp/(tp + fp)
  recall <- tp/(tp + fn) 
  f1 <- (2*precision*recall)/(precision + recall)
  
  return(f1)
}

## -----------------------------------------------------------------------------
tune_pick$F1 <- sapply(mxr_tune, function(x){
  compare <- merge(x, tune_ann, 
                   by = c("filename", "entity", "pos"), all = TRUE)
  prf(compare)
})


ggplot(tune_pick) + geom_point(aes(max_edit_distance, window_length, size = F1)) + 
  scale_y_continuous(breaks=seq(30,120,30)) + 
  annotate("text", x = tune_pick$max_edit_distance+.2, y = tune_pick$window_length,
           label = round(tune_pick$F1, 2)) + 
  ggtitle("F1 for tuning parameter values")

## ---- eval = FALSE------------------------------------------------------------
#  parseMedExtractR(filename)
#  parseMedXN(filename, begText = "^[ID0-9]+_[0-9-]+_[Note0-9]")
#  parseCLAMP(filename)
#  parseMedEx(filename)

## -----------------------------------------------------------------------------
tac_mxr_fn <- system.file("examples", "tac_mxr.csv", package = "EHR")
lam_mxr_fn <- system.file("examples", "lam_mxr.csv", package = "EHR")
lam_mxn_fn <- system.file("examples", "lam_medxn.csv", package = "EHR")
tac_mxr_parsed <- parseMedExtractR(tac_mxr_fn)
lam_mxr_parsed <- parseMedExtractR(lam_mxr_fn)

## ---- echo = FALSE, warning = FALSE-------------------------------------------
lam_mxn <- scan(lam_mxn_fn, '', sep = '\n', quiet = TRUE)
# Print output
message("MedXN output:\n")
#lam_mxn
cat(findbreaks(lam_mxn, charlen = 90))

## ---- warning = FALSE---------------------------------------------------------
lam_mxn_parsed <- parseMedXN(lam_mxn_fn, begText = "^[ID0-9]+_[0-9-]+_[Note0-9]")

## ---- echo = FALSE------------------------------------------------------------
cat(breaktable(lam_mxn_parsed, c(3,5)))

## ---- eval = FALSE------------------------------------------------------------
#  buildDose(dat, dn = NULL)

## ---- warning = FALSE---------------------------------------------------------
(tac_part_i_out <- buildDose(tac_mxr_parsed))
(lam_part_i_out <- buildDose(lam_mxr_parsed))

## -----------------------------------------------------------------------------
lam_checkForRare <- buildDose(lam_mxr_parsed, checkForRare=TRUE)

## -----------------------------------------------------------------------------
bmd <- function(x) {
  fns <- strsplit(x, '_')
  pid <- sapply(fns, `[`, 1)
  date <- as.Date(sapply(fns, `[`, 2), format = '%Y-%m-%d')
  note <- sapply(fns, `[`, 3)
  data.frame(filename = x, pid, date, note, stringsAsFactors = FALSE)
}
bmd("tacpid1_2008-06-26_note1_1.txt")
(tac_metadata <- bmd(tac_part_i_out[['filename']]))

## -----------------------------------------------------------------------------
data(tac_lab, package = 'EHR')
tac_lab

## -----------------------------------------------------------------------------
(tac_ld <- processLastDose(mxrData = tac_mxr, noteMetaData = tac_metadata, labData = tac_lab))

## -----------------------------------------------------------------------------
(tac_part_i_out_lastdose <- addLastDose(buildData = tac_part_i_out, lastdoseData = tac_ld))

## ---- eval = FALSE, warning = FALSE-------------------------------------------
#  collapseDose(x, noteMetaData, naFreq = 'most', ...)

## ---- eval = FALSE------------------------------------------------------------
#  tac_part_ii <- collapseDose(tac_part_i_out_lastdose, tac_metadata, naFreq = 'most')

## ---- echo = FALSE, warning = FALSE-------------------------------------------
suppressWarnings(tac_part_ii <- collapseDose(tac_part_i_out_lastdose, tac_metadata, naFreq = 'most'))

## -----------------------------------------------------------------------------
tac_part_ii$note

## -----------------------------------------------------------------------------
tac_part_ii$date

## ---- eval = FALSE------------------------------------------------------------
#  data(lam_metadata, package = 'EHR')
#  lam_part_ii <- collapseDose(lam_part_i_out, lam_metadata, naFreq = 'most', 'xr|er')

## ---- echo = FALSE, warning = FALSE-------------------------------------------
data(lam_metadata, package = 'EHR')
suppressWarnings(lam_part_ii <- collapseDose(lam_part_i_out, lam_metadata, naFreq = 'most', 'xr|er'))

## -----------------------------------------------------------------------------
lam_part_ii$note

## -----------------------------------------------------------------------------
lam_part_ii$date

## -----------------------------------------------------------------------------
x <- lam_part_ii[['note']]
# retrieve metadata for each filename
k1 <- lam_metadata[match(x[,'filename'], lam_metadata[,'filename']), c('pid','date','note')]
# select additional key data
k2 <- cbind(k1, x[,c('dose.daily','drugname_start')])
# turn keys into character string
chk <- do.call(paste, c(k2, sep = '|'))
# keep first instance of each chk key
lam_part_iii_note <- x[!duplicated(chk),]
lam_part_iii_note[,c('filename','drugname','drugname_start','dose.daily')]

## -----------------------------------------------------------------------------
x <- lam_part_ii[['date']]
# ignore note for date level collapsing
k1 <- lam_metadata[match(x[,'filename'], lam_metadata[,'filename']), c('pid','date')]
k2 <- cbind(k1, x[,c('dose.daily','drugname_start')])
chk <- do.call(paste, c(k2, sep = '|'))
lam_part_iii_date <- x[!duplicated(chk),]
lam_part_iii_date[,c('filename','drugname','drugname_start','dose.daily')]

