## Score the recovery studies by PSIS-LOO (and WAIC from the same pointwise values) exactly as by DIC, from the
## per-fit LOO results of loo.R (run jmh-loo1, n_draws 2000, DMC2 f2f58985).
##   CC_LOO_DIR=<folder of <spec>_<label>.rds> Rscript code/08_recovery/conflict/score_loo.R
## Inputs: recovery/conflict/recovery_jm.csv (studies 1-4), recovery_rec5.csv (study 5), recovery/choice/recch_dic.csv
##   (choice study 2), each with its fits' best-draw conditional DIC; <CC_LOO_DIR>/<spec>_<label>.rds (default
##   <ROOT>/recovery_work/loo). Read-only on all inputs.
## Outputs (write_out):
##   recovery/conflict/loo_recovery.csv  per fit: LOO, WAIC, DIC, p_loo, k counts, draws, minutes (all three studies)
##   recovery/conflict/loo_sets.csv      per data set x criterion: winner, recovered, gap = C(true) - C(winner),
##                                       se_gap (LOO/WAIC: 2 sqrt(n var(pointwise difference))), runner-up margin
##   recovery/conflict/loo_studies.csv   per study x criterion: sets, recovered, largest wrong-winner margin, and for
##                                       wrong winners whether their mechanism code is longer (more parameters) or shorter than the generator's
##   recovery/conflict/loo_rec5_cells.csv study 5 per cell x criterion: false alarms (G1 data won by a 2-mechanism
##                                       model), misses (G2 data won by G1) and their largest margins
##   recovery/choice/recch_loo.csv       choice sets: as loo_sets plus pseudo-BMA+ and stacking weights of the generator
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
LDIR <- Sys.getenv("CC_LOO_DIR", file.path(ROOT, "recovery_work", "loo"))
RC <- file.path(ROOT, "recovery", "conflict"); RCH <- file.path(ROOT, "recovery", "choice")
CRIT <- c("LOO", "WAIC", "DIC")
a <- read.csv(file.path(RC, "recovery_jm.csv"), stringsAsFactors = FALSE)
b <- read.csv(file.path(RC, "recovery_rec5.csv"), stringsAsFactors = FALSE)
ch <- read.csv(file.path(RCH, "recch_dic.csv"), stringsAsFactors = FALSE)
n_mech <- function(cand) nchar(gsub("[a-z]", "", sub("^[A-Z]+_", "", cand)))   # mechanisms = capital letters after the arch
cplx <- function(cand) nchar(sub("^[A-Z]+_", "", cand))   # complexity proxy: letters of the mechanism code (each adds parameters)
fits <- rbind(
  data.frame(study = a$study, set = sub("[.]RData$", "", a$data_file), spec = a$spec, candidate = a$candidate,
             generator = a$generator, true_model = a$true_model, DIC = a$DIC,
             label = ifelse(a$spec == "rXA_g2_r2_RDM_VciSci", paste0(a$spec, "_jmhrec1b"), paste0(a$spec, "_jmhrec1"))),
  data.frame(study = "S5", set = b$set, spec = b$spec, candidate = paste(b$arch, b$candidate, sep = "_"),
             generator = b$generator, true_model = b$true_model, DIC = b$DIC, label = paste0(b$spec, "_jmhrec5")),
  data.frame(study = "choice", set = sub("^(rc_[A-Z]+_r[0-9]+)_.*", "\\1", ch$spec), spec = ch$spec, candidate = ch$family,
             generator = sub("^Choice_([A-Z]+).*", "\\1", ch$generator), true_model = ch$true_model, DIC = ch$DIC,
             label = paste0(ch$spec, "_jmdrc2")))
fits$file <- file.path(LDIR, paste0(fits$label, ".rds")); have <- file.exists(fits$file)
cat(sum(have), "of", nrow(fits), "fits have LOO in", LDIR, "\n")
PW <- setNames(lapply(fits$file, function(f) if (file.exists(f)) readRDS(f) else NULL), fits$label)
S <- do.call(rbind, lapply(PW[have], `[[`, "summary"))
fits <- merge(fits, S[, c("label", "n_trials", "draws", "LOO", "WAIC", "p_loo", "p_waic", "se_elpd_loo", "mcse_elpd_loo",
                          "k_gt_0.7", "k_gt_1", "share_k_gt_0.7", "n_constant", "minutes", "emc2_sha")], by = "label", all.x = TRUE)
fits <- fits[order(fits$study, fits$set, fits$spec), ]
write_out(fits[, setdiff(names(fits), "file")], file.path(RC, "loo_recovery.csv"))

## per data set x criterion
pdiff <- function(l1, l2, what) { x <- PW[[l1]][[what]]; y <- PW[[l2]][[what]]
  if (is.null(x) || is.null(y) || length(x) != length(y)) return(NA_real_); 2 * sqrt(length(x) * var(x - y)) }
sets <- do.call(rbind, lapply(split(fits, fits$set), function(x) do.call(rbind, lapply(CRIT, function(cr) {
  v <- x[[cr]]; done <- !anyNA(v); o <- order(v); w <- o[1]; tr <- which(x$true_model)
  what <- c(LOO = "elpd_loo", WAIC = "elpd_waic", DIC = NA)[[cr]]
  data.frame(study = x$study[1], set = x$set[1], generator = x$generator[1], criterion = cr, n_cand = nrow(x), done = done,
             true = x$candidate[tr], winner = if (done) x$candidate[w] else NA, recovered = if (done) w == tr else NA,
             gap = if (done) v[tr] - v[w] else NA,
             se_gap = if (done && w != tr && !is.na(what)) pdiff(x$label[tr], x$label[w], what) else NA,
             margin = if (done) v[o[2]] - v[w] else NA,
             winner_cplx_minus_true = if (done && x$study[1] != "choice") cplx(x$candidate[w]) - cplx(x$candidate[tr]) else NA, row.names = NULL)
}))))
sets <- sets[order(sets$study, sets$set, sets$criterion), ]
write_out(sets, file.path(RC, "loo_sets.csv"))

studies <- do.call(rbind, lapply(split(sets, list(sets$study, sets$criterion), drop = TRUE), function(y) { y <- y[y$done, ]; if (!nrow(y)) return(NULL)
  wr <- y[!y$recovered, ]
  data.frame(study = y$study[1], criterion = y$criterion[1], sets = nrow(y), recovered = sum(y$recovered),
             max_wrong_margin = if (nrow(wr)) max(wr$gap) else 0,
             wrong_more_complex = sum(wr$winner_cplx_minus_true > 0), max_margin_more_complex = if (any(wr$winner_cplx_minus_true > 0)) max(wr$gap[wr$winner_cplx_minus_true > 0]) else 0,
             wrong_simpler = sum(wr$winner_cplx_minus_true < 0), max_margin_simpler = if (any(wr$winner_cplx_minus_true < 0)) max(wr$gap[wr$winner_cplx_minus_true < 0]) else 0,
             row.names = NULL) }))
studies <- studies[order(studies$study, match(studies$criterion, CRIT)), ]
write_out(studies, file.path(RC, "loo_studies.csv"))

## study 5 cells (G1 = one-mechanism generator, G2 = two-mechanism)
b$set <- sprintf("%s_%s_r%d", b$cell, b$gen, b$replicate)
cell_of <- unique(b[, c("set", "cell", "gen", "task", "arch")])
G1 <- tapply(b$generator[b$gen == "g1"], b$cell[b$gen == "g1"], `[`, 1)
s5 <- merge(sets[sets$study == "S5" & sets$done, ], cell_of, by = "set")
cells <- do.call(rbind, lapply(split(s5, list(s5$cell, s5$criterion), drop = TRUE), function(y) {
  g1 <- y[y$gen == "g1", ]; g2 <- y[y$gen == "g2", ]; g1name <- sub(paste0("^", y$task[1], "_"), "", G1[[y$cell[1]]])
  fa <- g1[n_mech(g1$winner) == 2, ]; mi <- g2[g2$winner == g1name, ]
  data.frame(cell = y$cell[1], task = y$task[1], arch = y$arch[1], criterion = y$criterion[1],
             sets_g1 = nrow(g1), recovered_g1 = sum(g1$recovered), false_alarms = nrow(fa), fa_max_margin = if (nrow(fa)) max(fa$gap) else 0,
             sets_g2 = nrow(g2), recovered_g2 = sum(g2$recovered), misses = nrow(mi), miss_max_margin = if (nrow(mi)) max(mi$gap) else 0,
             floor = max(y$gap), row.names = NULL) }))
cells <- cells[order(cells$cell, match(cells$criterion, CRIT)), ]
write_out(cells, file.path(RC, "loo_rec5_cells.csv"))

## choice sets: weights of the generator from the pointwise elpd
chs <- sets[sets$study == "choice", ]
wts <- do.call(rbind, lapply(split(fits[fits$study == "choice", ], fits$set[fits$study == "choice"]), function(x) {
  if (!all(x$label %in% names(PW)[have])) return(NULL)
  lp <- sapply(x$label, function(l) PW[[l]]$elpd_loo); set.seed(20261105)
  data.frame(set = x$set[1], w_pseudobma_true = as.numeric(loo::pseudobma_weights(lp, BB = TRUE))[x$true_model],
             w_stacking_true = as.numeric(loo::stacking_weights(lp))[x$true_model], row.names = NULL) }))
if (!is.null(wts)) chs <- merge(chs, wts, by = "set", all.x = TRUE)
write_out(chs, file.path(RCH, "recch_loo.csv"))

print(studies, row.names = FALSE, digits = 4)
cat("\nstudy 5 cells\n"); print(cells, row.names = FALSE, digits = 4)
cat("\nchoice\n"); print(chs[, intersect(names(chs), c("set", "criterion", "winner", "recovered", "gap", "se_gap", "w_pseudobma_true", "w_stacking_true"))], row.names = FALSE, digits = 4)
