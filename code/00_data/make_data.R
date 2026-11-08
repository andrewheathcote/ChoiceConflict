## Build the four analysed data sets from the published trial-level files of Eisenberg et al. (2019):
##   Rscript code/00_data/make_data.R [outdir]
##   reads inputs/raw/{choice_reaction_time,simon,stroop,attention_network_task}.csv
##   writes <outdir, default a temporary folder>/<task>/dat.RData (+ Stroop/dat_acc.RData) and compares each with the
##   archived inputs/data/<task>/ (all.equal on the data frames); pass outdir = inputs/data to overwrite them.
## Steps:
##   Choice  practice removed; participants with > 25% omissions, then omitted trials, then participants with accuracy
##           < .6 removed; trials with RT < 0.2 s removed (the original parseChoice.R).
##   Simon, Stroop, ANT   Data_Cleaning_Script.Rmd (participant exclusions: median RT <= 200 ms, > 25% omissions, accuracy
##           <= 60%, one key > 95%; trials with RT < 0.15 s or > mean + 3 SD removed; omissions rechecked). One fix, marked
##           in the script: its Stroop section reused the Simon section's excluded ids; the analysed Stroop data (which every
##           fit used) excludes only Stroop's own. The Rmd's other output (tables, plots, mixed models) is not used.
##   ANT     neutral flankers coded as congruent; cues centre + double -> alerting (the original ANT/all.R).
##   Stroop  dat_acc.RData: accuracy-coded responses (R = RS == R) for the diffusion models.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
a <- commandArgs(TRUE); OUT <- if (length(a)) a[1] else file.path(tempdir(), "data")
RAW <- file.path(INPUTS, "raw"); W <- file.path(tempdir(), "clean"); dir.create(file.path(W, "Data"), recursive = TRUE, showWarnings = FALSE)
for (f in c("choice_reaction_time", "simon", "stroop", "attention_network_task")) file.copy(file.path(RAW, paste0(f, ".csv")), file.path(W, "Data"))
save_dat <- function(dat, task, file = "dat.RData") { dir.create(file.path(OUT, task), recursive = TRUE, showWarnings = FALSE)
  save(dat, file = file.path(OUT, task, file)) }

## ---- simple choice ----------------------------------------------------------------------------------------------------
rawdat <- read.csv(file.path(RAW, "choice_reaction_time.csv"))
dat <- rawdat[rawdat$exp_stage != "practice", c(19, 4, 10, 3, 12)]
names(dat) <- c("subjects", "S", "R", "C", "rt"); dat <- dat[, -4]
dat$S <- factor(dat$S, labels = 1:2)
dat[dat$R == -1, c("R", "rt")] <- NA
dat$R <- factor(dat$R, labels = 1:2); dat$rt <- dat$rt / 1000
omit <- tapply(is.na(dat$R), dat$subjects, mean); dat <- dat[!(dat$subjects %in% names(omit)[omit > .25]), ]
dat <- dat[!is.na(dat$rt), ]
acc <- tapply(dat$S == dat$R, dat$subjects, mean); dat <- dat[!(dat$subjects %in% names(acc)[acc < .6]), ]
dat$subjects <- factor(dat$subjects); dat <- dat[!(dat$rt < .2), ]
save_dat(dat, "Choice")

## ---- conflict tasks: the cleaning script ------------------------------------------------------------------------------
knitr::purl(file.path(CODE, "00_data", "Data_Cleaning_Script.Rmd"), output = file.path(W, "clean.R"), quiet = TRUE)
owd <- setwd(W); ev <- new.env(); sys.source(file.path(W, "clean.R"), envir = ev); setwd(owd)
for (t in c("simon", "stroop")) {
  d <- as.data.frame(get(paste0("data_", t), ev)); rownames(d) <- NULL
  save_dat(d, if (t == "simon") "Simon" else "Stroop")
}
s <- as.data.frame(get("data_stroop", ev)); rownames(s) <- NULL
s$R <- factor(s$RS == s$R); dat <- s; save_dat(dat, "Stroop", "dat_acc.RData")

## ---- ANT recode -------------------------------------------------------------------------------------------------------
ant <- as.data.frame(get("data_ant", ev))[, c("subject", "rt", "stim_relevant", "stim_irrelevant", "congruence", "cue", "R")]
ant$stim_flanker <- ifelse(ant$congruence == "congruent", ant$stim_relevant,
                    ifelse(ant$congruence == "incongruent", ifelse(ant$stim_relevant == "left", "right", "left"), "neutral"))
ant <- ant[, c("subject", "rt", "stim_relevant", "stim_flanker", "stim_irrelevant", "cue", "R")]
colnames(ant) <- c("subjects", "rt", "RS", "IS", "IS_loc", "CUE", "R")
for (k in c("subjects", "R", "RS", "IS", "IS_loc", "CUE")) ant[[k]] <- as.factor(ant[[k]])
dat <- ant
dat$IS <- as.character(dat$IS); dat$IS[dat$IS == "neutral"] <- as.character(dat$RS[dat$IS == "neutral"]); dat$IS <- factor(dat$IS)
names(dat)[5] <- "L"
levels(dat$CUE) <- c("alert", "alert", "nocue", "orient"); names(dat)[6] <- "Q"
dat$Q <- factor(as.character(dat$Q), levels = c("nocue", "alert", "orient"))
save_dat(dat, "ANT")

## ---- compare with the archived data -------------------------------------------------------------------------------------
for (x in list(c("Choice", "dat.RData"), c("Simon", "dat.RData"), c("Stroop", "dat.RData"), c("Stroop", "dat_acc.RData"), c("ANT", "dat.RData"))) {
  new <- local({ e <- new.env(); load(file.path(OUT, x[1], x[2]), envir = e); e$dat })
  old <- load_dat(x[1], x[2])
  cat(sprintf("%-7s %-14s rows %6d / %6d  participants %3d / %3d  equal: %s\n", x[1], x[2], nrow(new), nrow(old),
              nlevels(factor(new$subjects)), nlevels(factor(old$subjects)), isTRUE(all.equal(new, old, check.attributes = FALSE))))
}
