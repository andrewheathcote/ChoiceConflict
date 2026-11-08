## ChoiceConflict: paths and EMC2 builds, shared by every script. Nothing here is machine-specific.
## A script in code/<stage>/ starts with
##   CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
##   source(file.path(CODE, "config.R"))
## Environment (all optional):
##   CC_ROOT     the archive root holding inputs/, fits/, derived/ (default: the folder above code/)
##   CC_LIB      where install_builds.sh put the builds, one library per sha: <CC_LIB>/<sha8>/EMC2 (default <CC_ROOT>/lib)
##   CC_FIGURES  where figures are written (default <CC_ROOT>/figures)
##   CC_DERIVED  where derived tables and predictives are read and written (default <CC_ROOT>/derived)
##   CC_INPUTS   specs and data (default <CC_ROOT>/inputs; the recovery studies use <work>/inputs, see 08_recovery)
##   CC_REFITS   where 02_fit/fit.R writes (default <CC_ROOT>/refits)
if (!exists("CODE")) stop("set CODE (the code/ folder) before sourcing config.R")
ROOT    <- normalizePath(Sys.getenv("CC_ROOT", dirname(CODE)), mustWork = TRUE)
LIBDIR  <- Sys.getenv("CC_LIB", file.path(ROOT, "lib"))
FIGDIR  <- Sys.getenv("CC_FIGURES", file.path(ROOT, "figures"))
DERIVED <- Sys.getenv("CC_DERIVED", file.path(ROOT, "derived"))
INPUTS  <- Sys.getenv("CC_INPUTS", file.path(ROOT, "inputs"))   # recovery studies point this at their own specs and data
TABLES  <- file.path(DERIVED, "tables")
TASKS   <- c("Simon", "Stroop", "ANT")

## the archive's own tables and folders
tab_path  <- function(name) file.path(TABLES, name)
read_tab  <- function(name) read.csv(tab_path(name), stringsAsFactors = FALSE)
spec_path <- function(spec) file.path(INPUTS, "specs", paste0(spec, ".rds"))
read_spec <- function(spec) readRDS(spec_path(spec))
## data: inputs/data/<task>/dat.RData; Stroop's diffusion models use the accuracy-coded dat_acc.RData
dat_path  <- function(task, file = "dat.RData") file.path(INPUTS, "data", task, file)
load_dat  <- function(task, file = "dat.RData") { ev <- new.env(); load(dat_path(task, file), envir = ev); ev$dat }
spec_dat  <- function(sp) load_dat(sp$task, basename(sp$data_file))
## fits: fits/<task>/<arch>/<label>.RData (FITS.csv: label, spec, task, arch, build, class)
fit_path  <- function(label) {
  f <- Sys.glob(file.path(ROOT, "fits", "*", "*", paste0(label, ".RData")))
  if (length(f) != 1) stop("fit ", label, " not in ", file.path(ROOT, "fits"), " (", length(f), " matches)")
  f
}
load_fit  <- function(label) { ev <- new.env(); nm <- load(fit_path(label), envir = ev); ev[[nm[1]]] }
pp_subj_path <- function(label) file.path(DERIVED, "pp_subj", paste0(label, ".rds"))
pp_raw_path  <- function(label) file.path(DERIVED, "pp_raw", paste0(label, "_subj.rds"))
fig_path  <- function(name) { dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE); file.path(FIGDIR, name) }

## EMC2 builds (builds/BUILDS.csv). Each stage loads the build that made its objects; use_build() stops if the
## library is missing or is a different commit.
use_build <- function(sha) {
  lib <- file.path(LIBDIR, substr(sha, 1, 8))
  if (!dir.exists(file.path(lib, "EMC2"))) stop("EMC2 build ", sha, " is not installed in ", lib, ": run builds/install_builds.sh")
  got <- readLines(file.path(lib, "EMC2_SHA.txt"))
  if (!startsWith(got, substr(sha, 1, 8))) stop("library ", lib, " holds EMC2 ", got, ", not ", sha)
  suppressMessages(library(EMC2, lib.loc = lib))
  invisible(got)
}
## the builds by role (see builds/BUILDS.csv for which fits each made)
BUILD <- c(predict = "3c1bcaf1",    # subject-level predictives (pp_subj, pp_raw) for the conflict fits
           choice  = "4f1fa123",    # simple-choice fits, their DIC and predictives
           ic      = "f2f58985",    # WAIC / PSIS-LOO / pseudo-BMA+ / stacking (choice and its recovery study)
           jm      = "d0b6cf0b")    # the joint-move refits (label jm) and the Q6 restricted Simon LNR refits (jmr)

## quiet evaluation of chatty package calls
quietly <- function(x) { f <- tempfile(); sink(f); on.exit({ sink(); unlink(f) }); force(x) }
## write a csv only inside the archive root (never elsewhere)
write_out <- function(x, file) {
  stopifnot(startsWith(normalizePath(dirname(file), mustWork = FALSE), ROOT) || startsWith(file, DERIVED) || startsWith(file, FIGDIR))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file, row.names = FALSE); invisible(file)
}
