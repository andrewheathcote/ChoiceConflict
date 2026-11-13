## Population-mean conflict effects, one row per model (RDM_ViSi, RDM_BiSi, LNR_MiSi), one bar per task: bar height the
## posterior median, whisker the 95% credible interval (2.5% and 97.5%), from the pooled sample-stage mu draws of the
## displayed fit. Four panels per row (EEA, St, then Vi / Bi / Mi, then Si), each effect a difference match - mismatch
## (Andrew, 2026-11-10/11). One y scale, across the whole figure, per unit and kind: EEA with the i mechanism (Vi, Mi;
## Bi, in E, is alone) and St with Si:
##   EEA  evidence-accumulation advantage: accumulator matching the relevant attribute minus the other (lM); RDM v, LNR m
##   i    the irrelevant-match mechanism: accumulator matching the irrelevant attribute minus the other (iM);
##        RDM Vi v, RDM Bi B, LNR Mi m
##   St   variability, lM: RDM s, LNR s
##   Si   variability, iM: RDM s, LNR s
## Each effect is a difference on the natural scale, with every other contrast at zero (the design's grand mean):
## RDM v, B, s and LNR s are sampled on the log scale (RDM s intercept fixed at log 1), LNR m on the natural scale, and
## v_iMd, B_iMd, m_iMd are order restricted (exp). Contrast signs are read from each spec's clist, not assumed.
## Units: RDM rates and their variability E/s (E = evidence unit, s = seconds), RDM thresholds E, LNR log(s).
##   Rscript code/07_figures/plot_effects.R [MODEL ...]    (default RDM_ViSi RDM_BiSi LNR_MiSi)
## Reads derived/tables/display_jm.csv, inputs/specs/<spec>.rds and the fits; writes derived/tables/effects_mu.csv
## (plot-ready) and figures/effects.pdf.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
MODELS <- commandArgs(trailingOnly = TRUE); if (!length(MODELS)) MODELS <- c("RDM_ViSi", "RDM_BiSi", "LNR_MiSi")
D <- read_tab("display_jm.csv")

mu_draws <- function(label) { f <- load_fit(label)
  do.call(rbind, lapply(f, function(ch) t(ch$samples$theta_mu[, ch$samples$stage == "sample", drop = FALSE]))) }
## contrast weight of level `lev` of the logical factor `fac` (levels FALSE, TRUE: the clist rows in that order)
w <- function(sp, p, fac, lev) { m <- sp$clist[[p]][[fac]]; stopifnot(dim(m) == c(2, 1)); m[match(lev, c("FALSE", "TRUE")), 1] }
## natural-scale difference x(match) - x(mismatch) for parameter p, effect column `col` (exp'd if restricted)
diff_nat <- function(M, sp, p, col, fac, log_scale, int = M[, p]) {
  b <- M[, col]; if (col %in% names(sp$pre$func)) b <- exp(b)
  eta <- function(lev) int + w(sp, p, fac, lev) * b
  if (log_scale) exp(eta("TRUE")) - exp(eta("FALSE")) else eta("TRUE") - eta("FALSE") }

effects <- function(spec) {
  sp <- readRDS(file.path(ROOT, "inputs", "specs", paste0(spec, ".rds"))); M <- mu_draws(D$label[D$spec == spec])
  rdm <- sp$model == "RDM"; r <- if (rdm) "v" else "m"
  i <- sub("^.*_([A-Z]i)Si$", "\\1", spec); ip <- c(Vi = "v", Bi = "B", Mi = "m")[[i]]
  s_int <- if (rdm) { stopifnot(sp$consts[["s"]] == 0); rep(0, nrow(M)) } else M[, "s"]   # RDM s intercept fixed at log 1
  e <- list(EEA = diff_nat(M, sp, r, paste0(r, "_lMd"), "lM", rdm),
            I   = diff_nat(M, sp, ip, paste0(ip, "_iMd"), "iM", rdm),
            St  = diff_nat(M, sp, "s", "s_lMd", "lM", TRUE, int = s_int),
            Si  = diff_nat(M, sp, "s", "s_iMd", "iM", TRUE, int = s_int))
  names(e)[2] <- i
  cat(sprintf("%-16s %s  %d draws\n", spec, D$label[D$spec == spec], nrow(M)))
  sapply(e, quantile, c(0.025, 0.5, 0.975)) }

unit_of <- function(model, k) if (grepl("^LNR", model)) "log(s)" else if (k == "Bi") "E" else "E/s"
T <- do.call(rbind, lapply(MODELS, function(MODEL) do.call(rbind, lapply(TASKS, function(TK) {
  x <- effects(paste(TK, MODEL, sep = "_"))
  data.frame(model = MODEL, task = TK, effect = colnames(x), unit = sapply(colnames(x), unit_of, model = MODEL),
             lo = x[1, ], median = x[2, ], hi = x[3, ], row.names = NULL) }))))
write_out(T, tab_path("effects_mu.csv"))

T$scale <- paste(T$unit, ifelse(T$effect %in% c("St", "Si"), "variability", "location"))
pad <- function(r) r + c(-0.04, 0.04) * diff(r)
YL <- lapply(split(T, T$scale), function(x) pad(range(0, x$lo, x$hi)))
## fixed ranges (Andrew, 2026-11-11), by "<arch> <effect>"; EEA and St keep YL
FIXED <- list("RDM Vi" = c(0, 1), "RDM Bi" = c(-1, 0), "RDM Si" = c(-0.25, 0.25), "LNR Mi" = c(-0.25, 0), "LNR Si" = c(-0.15, 0.15))
FONT <- "Palatino"                                       # model and mechanism names: an i that cannot pass for an l
OUT <- fig_path("effects.pdf")
pdf(OUT, width = 12, height = 3.2 * length(MODELS))
par(mfrow = c(length(MODELS), 4), mar = c(2.4, 4.2, 2, 0.6), oma = c(0, 2.6, 0, 0), mgp = c(2.6, 0.5, 0), tcl = -0.25, las = 1)
for (MODEL in MODELS) { x0 <- T[T$model == MODEL, ]; arch <- sub("_.*", "", MODEL)
  for (k in c("EEA", "St", setdiff(unique(x0$effect), c("EEA", "St", "Si")), "Si")) {
    Q <- x0[x0$effect == k, ]; Q <- Q[match(TASKS, Q$task), ]
    x <- barplot(Q$median, ylim = if (paste(arch, k) %in% names(FIXED)) pad(FIXED[[paste(arch, k)]]) else YL[[Q$scale[1]]], col = "grey80", border = "grey30", names.arg = TASKS, ylab = Q$unit[1],
                 cex.names = 1.1, cex.axis = 0.9, cex.lab = 1.1, width = 1, space = 0.4)
    abline(h = 0, col = "grey30")
    arrows(x, Q$lo, x, Q$hi, angle = 90, code = 3, length = 0.05, lwd = 1.2)
    title(main = k, cex.main = 1.4, line = 0.6, family = FONT)
    box(bty = "l")
    if (k == "EEA") mtext(sub("_", " ", MODEL), side = 2, line = 4.9, cex = 1.15, font = 2, las = 0, family = FONT) } }
invisible(dev.off())
if (nzchar(Sys.which("gs"))) embedFonts(OUT) else warning("Ghostscript not found: fonts not embedded")   # Palatino is not a base PDF font
cat("wrote", OUT, "\n")
print(transform(T, lo = signif(lo, 3), median = signif(median, 3), hi = signif(hi, 3)), row.names = FALSE)
