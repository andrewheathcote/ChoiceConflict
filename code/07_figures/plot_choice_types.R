## Simple choice: delta and CAF by stimulus for the DIC-best variant of each architecture other than the DDM (Figure 4 /
## choice.pdf draws the DDM's), one row per architecture (WDM, RDM, LNR, LBA), from the same subject-level predictives
## as choice.pdf (derived/pp_raw/<label>_subj.rds, made by plot_choice.R sim) and EMC2's plot_delta / plot_caf with
## choice.pdf's settings.
##   Rscript code/07_figures/plot_choice_types.R          figures/choice_types.pdf
## Inputs (read-only): inputs/data/Choice/, derived/tables/choice_dic.csv, derived/pp_raw/.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
use_build(BUILD[["choice"]])
OUT <- fig_path("choice_types.pdf")
TYPES <- c("WDM", "RDM", "LNR", "LBA")
## each architecture's DIC-best variant; its predictives exist only for the full variants (CHOICE_FULL)
D <- read_tab("choice_dic.csv"); D <- D[order(D$DIC), ]
fam <- sub("^(LBA|RDM|LNR|WDM|DDM).*$", "\\1", D$variant)
best <- setNames(D$label[match(TYPES, fam)], TYPES)
stopifnot(identical(unname(best), unname(CHOICE_FULL[TYPES])))
dat <- load_dat("Choice")
pp <- lapply(best, function(lb) readRDS(pp_raw_path(lb)))
pct <- c(10, 50, 90)
ttl <- c(WDM = "Wiener Diffusion Model (WDM)", RDM = "Racing Diffusion Model (RDM)", LNR = "Lognormal Race (LNR)",
         LBA = "Linear Ballistic Accumulator (LBA)")
pdf(OUT, width = 8, height = 14)
par(mfrow = c(4, 2))
for (m in TYPES) {
  plot_delta(dat, pp[[m]], delta_factor = "Stimulus", main = paste(ttl[[m]], "\nDelta Function"), layout = NULL, legendpos = NA,
             add_percentiles = pct, functions = list(Stimulus = \(d)d$S), ylab = "RT (2-1) (sec)")
  plot_caf(dat, pp[[m]], caf_factor = "Stimulus", main = paste(ttl[[m]], "\nConditional Accuracy Function"), layout = NULL,
           ylab = "Accuracy (%)", smooth_window = 10, legendpos = c("bottomright", NA), functions = list(Stimulus = \(d)d$S),
           accuracy_function = function(d) d$S == d$R, ylim = c(90, 100))
}
invisible(dev.off())
cat("wrote", OUT, "\n")
