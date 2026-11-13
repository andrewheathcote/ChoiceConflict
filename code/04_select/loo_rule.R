## A LOO version of the paper's DIC "=" rule, calibrated on the recovery studies the way the DIC thresholds were, and
## applied to the real fits, to say how model selection would differ under LOO. Read-only on all inputs.
##   CC_LOO_DIR=<recovery pointwise folder> Rscript code/04_select/loo_rule.R
## DIC rule (plot_heat_grid.R, threshold_compare.R): within an architecture a model is a contender ("=") with the best
## if it trails by <= ADD (50) when it has FEWER mechanisms than the best, else by <= DROP (150). ADD and DROP are the
## largest margins by which a wrong winner beat the generator in recovery studies 1-5 (more mechanisms than the
## generator: 40.4, 27.8; fewer or as many: 149.4), rounded up to the next 10. Mechanisms = tokens (Vci = Vc + Vi).
## LOO rules, calibrated identically on LOO's recovery margins (studies 1-5 incl. cross-architecture):
##   LOO_dir    directional ADD_L / DROP_L                       (the analogue of the paper rule)
##   LOO_single one floor DROP_L both ways                         (the analogue of the single 150 floor)
##   LOO_2SE    contender if it trails by <= 2 SE of the pointwise difference (usual LOO practice; not calibrated)
##   LOO_dir2SE contender if within LOO_dir OR within 2 SE
## A mechanism is "needed" (threshold_compare.R section 2) if the best model without it trails the best by more than
## the rule's threshold for that comparison.
## Inputs: recovery/conflict/loo_recovery.csv (score_loo.R), <CC_LOO_DIR>/<label>.rds, derived/tables/loo_real.csv,
##   derived/loo_real/<label>.rds.
## Outputs: recovery/conflict/loo_rule_calib.csv   thresholds per criterion (raw max margins, rounded; second largest), with variants
##          recovery/conflict/loo_rule_recovery.csv per study x rule: coverage of the generator, decisive sets,
##                                                   contenders, parsimonious pick; study 5 false alarms / detections
##          derived/tables/loo_rule_within.csv      per task x architecture x rule: best, contenders, mechanisms needed
##          derived/tables/loo_rule_mechs.csv       per task x architecture x rule x mechanism of the best: needed?
##          derived/tables/loo_rule_arch.csv        per task x rule: architectures not separable from the best one
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
LDIR <- Sys.getenv("CC_LOO_DIR", file.path(ROOT, "recovery_work", "loo")); RDIR <- file.path(DERIVED, "loo_real")
RC <- file.path(ROOT, "recovery", "conflict")
ADD_DIC <- 50; DROP_DIC <- 150
up10 <- function(x) ceiling(x / 10) * 10
tok <- function(v) { m <- regmatches(v, gregexpr("(SV|[A-Z])[ci]*", v))[[1]]
  unlist(lapply(m, function(t) { p <- sub("[ci]+$", "", t); s <- sub("^(SV|[A-Z])", "", t); paste0(p, strsplit(s, "")[[1]]) })) }
PWC <- new.env()
pw <- function(dir, l) { k <- paste(dir, l); if (is.null(PWC[[k]])) PWC[[k]] <- readRDS(file.path(dir, paste0(l, ".rds")))$elpd_loo; PWC[[k]] }
se_d <- function(dir, l1, l2) {   # SE of a LOO (deviance-scale) difference
 if (l1 == l2) return(0); x <- pw(dir, l1); y <- pw(dir, l2); 2 * sqrt(length(x) * var(x - y)) }

## contenders of the best under a rule: x has var, n (mechanisms), crit value v, label; returns logical
RULES <- c("DIC_dir", "DIC_single", "LOO_dir", "LOO_single", "LOO_2SE", "LOO_dir2SE")
crit_of <- function(rule) sub("_.*", "", rule)
contend <- function(x, rule, th, dir) {
  v <- x[[crit_of(rule)]]; b <- which.min(v); g <- v - v[b]; fewer <- x$n < x$n[b]
  se <- if (grepl("2SE", rule)) 2 * vapply(x$label, function(l) se_d(dir, x$label[b], l), 0) else NA   # 2 SE
  thr <- switch(rule, DIC_dir = ifelse(fewer, ADD_DIC, DROP_DIC), DIC_single = DROP_DIC,
                LOO_dir = ifelse(fewer, th$ADD, th$DROP), LOO_single = th$DROP, LOO_2SE = se,
                LOO_dir2SE = pmax(ifelse(fewer, th$ADD, th$DROP), se))
  list(best = b, gap = g, thr = thr, cont = g <= thr)
}

## ---------------------------------------------------------------- recovery: calibrate
F <- read.csv(file.path(RC, "loo_recovery.csv"), stringsAsFactors = FALSE); F <- F[F$study != "choice", ]
F$var <- sub("^[A-Z]+_", "", F$candidate); F$n <- lengths(lapply(F$var, tok))
F$rdm_old <- F$study == "S5" & grepl("_RDM_", F$label)          # study 5 RDM cells: sampled with the old RDM likelihood
S <- split(F, F$set)
margins <- do.call(rbind, lapply(S, function(x) do.call(rbind, lapply(c("LOO", "DIC"), function(cr) {
  b <- which.min(x[[cr]]); t <- which(x$true_model)
  data.frame(study = x$study[1], set = x$set[1], criterion = cr, gap = x[[cr]][t] - x[[cr]][b],
             truth_fewer = x$n[t] < x$n[b], rdm_old = x$rdm_old[1], row.names = NULL) }))))
cal <- function(m, what) do.call(rbind, lapply(c("LOO", "DIC"), function(cr) { y <- m[m$criterion == cr, ]
  a <- max(c(0, y$gap[y$truth_fewer])); d <- max(c(0, y$gap[!y$truth_fewer]))
  s2 <- function(g) if (length(g) > 1) sort(g, decreasing = TRUE)[2] else 0
  data.frame(criterion = cr, variant = what, max_gap_complex_wins = a, max_gap_other = d, ADD = up10(a), DROP = up10(d),
             second_gap_complex_wins = s2(y$gap[y$truth_fewer]), second_gap_other = s2(y$gap[!y$truth_fewer])) }))
calib <- rbind(cal(margins, "studies 1-5"), cal(margins[margins$study != "XA", ], "without cross-architecture"),
               cal(margins[!margins$rdm_old, ], "without study 5 RDM cells"),
               cal(margins[margins$study != "S5", ], "studies 1-4 only"))
write_out(calib, file.path(RC, "loo_rule_calib.csv"))
TH <- list(ADD = calib$ADD[calib$criterion == "LOO" & calib$variant == "studies 1-5"],
           DROP = calib$DROP[calib$criterion == "LOO" & calib$variant == "studies 1-5"])
stopifnot(calib$ADD[calib$criterion == "DIC" & calib$variant == "studies 1-5"] == ADD_DIC,
          calib$DROP[calib$criterion == "DIC" & calib$variant == "studies 1-5"] == DROP_DIC)   # reproduces the paper's 50 / 150
cat(sprintf("LOO thresholds: ADD %g, DROP %g (DIC: %g, %g)\n", TH$ADD, TH$DROP, ADD_DIC, DROP_DIC)); print(calib, row.names = FALSE)

## ---------------------------------------------------------------- recovery: apply each rule
RS <- do.call(rbind, lapply(S, function(x) do.call(rbind, lapply(RULES, function(r) {
  k <- contend(x, r, TH, LDIR); t <- which(x$true_model); cn <- which(k$cont)
  pick <- cn[order(x$n[cn], x[[crit_of(r)]][cn])][1]                 # parsimonious pick: fewest mechanisms, then best value
  one <- x$n == 1
  data.frame(study = x$study[1], set = x$set[1], gen = if (x$study[1] == "S5") sub(".*_(g[12])_r[0-9]+$", "\\1", x$set[1]) else NA,
             rule = r, n_cand = nrow(x), covered = k$cont[t], decisive = length(cn) == 1 && cn == t, n_cont = length(cn),
             pick_ok = pick == t, best_ok = k$best == t,
             second_needed = if (x$study[1] == "S5") !any(k$cont & one) else NA, row.names = NULL) }))))
RS <- RS[order(RS$study, RS$set, match(RS$rule, RULES)), ]
rec <- do.call(rbind, lapply(split(RS, list(RS$study, RS$rule), drop = TRUE), function(y) data.frame(
  study = y$study[1], rule = y$rule[1], sets = nrow(y), best_recovered = sum(y$best_ok), covered = sum(y$covered),
  decisive = sum(y$decisive), mean_contenders = mean(y$n_cont), pick_recovered = sum(y$pick_ok),
  s5_false_alarms = if (y$study[1] == "S5") sum(y$second_needed[y$gen == "g1"]) else NA,
  s5_detected = if (y$study[1] == "S5") sum(y$second_needed[y$gen == "g2"]) else NA, row.names = NULL)))
tot <- do.call(rbind, lapply(split(RS, RS$rule), function(y) data.frame(study = "all", rule = y$rule[1], sets = nrow(y),
  best_recovered = sum(y$best_ok), covered = sum(y$covered), decisive = sum(y$decisive), mean_contenders = mean(y$n_cont),
  pick_recovered = sum(y$pick_ok), s5_false_alarms = NA, s5_detected = NA)))
rec <- rbind(rec, tot); rec <- rec[order(match(rec$study, c("LBA", "LNR", "RDM", "XA", "S5", "all")), match(rec$rule, RULES)), ]
write_out(rec, file.path(RC, "loo_rule_recovery.csv"))
cat("\nrecovery under each rule\n"); print(rec, row.names = FALSE, digits = 3)

## ---------------------------------------------------------------- real fits
R <- read_tab("loo_real.csv"); R <- R[!is.na(R$DIC) & R$display_converged, ]
R$var <- sub("^[^_]+_[^_]+_", "", R$spec); R$n <- lengths(lapply(R$var, tok))
W <- NULL; M <- NULL
for (x in split(R, list(R$task, R$arch), drop = TRUE)) for (r in RULES) {
  k <- contend(x, r, TH, RDIR); b <- k$best; cn <- which(k$cont); v <- x[[crit_of(r)]]
  le2 <- which(x$n <= 2); b2 <- le2[which.min(v[le2])]
  pick <- cn[order(x$n[cn], v[cn])][1]
  ## mechanisms of the best: needed if the best model without it is not a contender
  need <- character(0)
  for (m in tok(x$var[b])) { wo <- which(!vapply(lapply(x$var, tok), function(t) m %in% t, TRUE)); if (!length(wo)) next
    w <- wo[which.min(v[wo])]
    M <- rbind(M, data.frame(task = x$task[1], arch = x$arch[1], rule = r, best = x$var[b], mech = m, best_without = x$var[w],
                             gap = v[w] - v[b], threshold = k$thr[w], needed = !k$cont[w], row.names = NULL))
    if (!k$cont[w]) need <- c(need, m) }
  W <- rbind(W, data.frame(task = x$task[1], arch = x$arch[1], rule = r, shown = nrow(x), best = x$var[b], best_n = x$n[b],
             n_contenders = length(cn), contenders = paste(sprintf("%s+%.0f", x$var[setdiff(cn, b)], k$gap[setdiff(cn, b)]), collapse = " "),
             min_mech_contender = min(x$n[cn]), pick = x$var[pick], best_le2 = if (length(b2)) x$var[b2] else NA,
             best_le2_gap = if (length(b2)) k$gap[b2] else NA,
             best_le2_gap_per_trial = if (length(b2)) k$gap[b2] / (2 * x$n_trials[1]) else NA,   # nats per trial
             two_ok = length(b2) > 0 && k$cont[b2],
             needed = paste(need, collapse = " "), n_needed = length(need), row.names = NULL)) }
W <- W[order(W$task, W$arch, match(W$rule, RULES)), ]
write_out(W, tab_path("loo_rule_within.csv")); write_out(M, tab_path("loo_rule_mechs.csv"))
cat("\nreal fits within architecture\n")
print(W[, c("task", "arch", "rule", "best", "n_contenders", "min_mech_contender", "pick", "best_le2", "best_le2_gap", "two_ok", "needed")], row.names = FALSE, digits = 4)

## architectures: best model of each; others are contenders of the best architecture within DROP (or 2 SE)
A <- NULL
for (x in split(R, R$task)) for (r in RULES) { cr <- crit_of(r)
  bb <- do.call(rbind, lapply(split(x, x$arch), function(y) y[which.min(y[[cr]]), ])); bb <- bb[order(bb[[cr]]), ]
  g <- bb[[cr]] - bb[[cr]][1]
  thr <- switch(r, DIC_dir = , DIC_single = rep(DROP_DIC, nrow(bb)), LOO_dir = , LOO_single = rep(TH$DROP, nrow(bb)),
                LOO_2SE = 2 * vapply(bb$label, function(l) se_d(RDIR, bb$label[1], l), 0),
                LOO_dir2SE = pmax(TH$DROP, 2 * vapply(bb$label, function(l) se_d(RDIR, bb$label[1], l), 0)))
  gn <- c(diff(bb[[cr]]), NA)
  sen <- if (grepl("2SE", r)) 2 * c(vapply(seq_len(nrow(bb) - 1), function(i) se_d(RDIR, bb$label[i], bb$label[i + 1]), 0), NA) else NA
  A <- rbind(A, data.frame(task = x$task[1], rule = r, order = paste(bb$arch, collapse = " < "),
             gaps = paste(sprintf("%s %.0f", bb$arch, g), collapse = ", "),
             contenders_of_best = paste(bb$arch[g <= thr & g > 0], collapse = " "),
             unresolved_adjacent = paste(na.omit(ifelse(gn <= switch(r, DIC_dir = , DIC_single = DROP_DIC, LOO_dir = , LOO_single = TH$DROP,
                                                                    LOO_2SE = sen, LOO_dir2SE = pmax(TH$DROP, sen)),
                                                       paste0(bb$arch, "~", c(bb$arch[-1], NA)), NA)), collapse = " "), row.names = NULL)) }
write_out(A, tab_path("loo_rule_arch.csv"))
cat("\narchitectures\n"); print(A, row.names = FALSE)

## reference: in study 5 two-mechanism (g2) data, how much the generating second mechanism gains over the best
## one-mechanism model, by each criterion (deviance scale; per trial in nats)
g2 <- do.call(rbind, lapply(S[grepl("_g2_", names(S))], function(x) { t <- which(x$true_model); o <- which(x$n == 1); if (!length(o)) return(NULL)
  data.frame(set = x$set[1], LOO = min(x$LOO[o]) - x$LOO[t], DIC = min(x$DIC[o]) - x$DIC[t], per_trial_LOO = (min(x$LOO[o]) - x$LOO[t]) / (2 * x$n_trials[1])) }))
cat("\nstudy 5 g2 data: gain of the generator over the best one-mechanism model\n"); print(sapply(g2[, -1], quantile, c(0, .1, .25, .5, .75, 1)), digits = 3); cat(nrow(g2), "sets\n")
