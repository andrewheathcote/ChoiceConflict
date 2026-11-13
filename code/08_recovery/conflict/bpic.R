## Why the conflict tasks are scored by DIC and not BPIC: every conflict recovery study rescored by BPIC (Ando, 2007).
##   Rscript code/08_recovery/conflict/bpic.R    (reads recovery/conflict/recovery_jm.csv and recovery_rec5.csv)
## BPIC = meanD + 2 pD with the best-draw pD = DIC - meanD (the per-subject minimum is already in the stored DIC), so
## BPIC = 2 DIC - meanD. Per data set and criterion: the winner, whether it is the generator, and the true-vs-winner gap.
## Per study (and, for study 5, per generator: G1 one mechanism, G2 two): sets recovered and the largest wrong-winner gap.
## Writes recovery/conflict/bpic_sets.csv and bpic_recovery.csv.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
RC <- file.path(ROOT, "recovery", "conflict")
a <- read.csv(file.path(RC, "recovery_jm.csv"), stringsAsFactors = FALSE)
a$set <- sprintf("%s_%s_r%d", a$study, a$gen, a$replicate); a$group <- paste("study", match(a$study, c("LNR", "RDM", "LBA", "XA")))
b <- read.csv(file.path(RC, "recovery_rec5.csv"), stringsAsFactors = FALSE)
b$group <- ifelse(b$gen == "g1", "study 5, G1 data", "study 5, G2 data")
keep <- c("group", "set", "generator", "candidate", "true_model", "DIC", "meanD")
m <- rbind(a[, keep], b[, keep]); stopifnot(!anyNA(m$DIC), !anyNA(m$meanD))
m$BPIC <- 2 * m$DIC - m$meanD
S <- do.call(rbind, lapply(split(m, m$set), function(x) do.call(rbind, lapply(c("DIC", "BPIC"), function(cr) {
  d <- setNames(x[[cr]], x$candidate); w <- which.min(d)
  data.frame(criterion = cr, group = x$group[1], set = x$set[1], generator = x$generator[1], winner = names(d)[w],
             recovered = x$true_model[w], gap = round(d[x$true_model] - min(d), 1), row.names = NULL) }))))
S <- S[order(S$criterion != "DIC", S$group, S$set), ]
write_out(S, file.path(RC, "bpic_sets.csv"))
R <- do.call(rbind, lapply(split(S, list(S$criterion, S$group), drop = TRUE), function(s)
  data.frame(criterion = s$criterion[1], group = s$group[1], sets = nrow(s), recovered = sum(s$recovered),
             largest_wrong_winner_gap = max(s$gap))))
R <- R[order(R$group, R$criterion != "DIC"), ]; rownames(R) <- NULL
write_out(R, file.path(RC, "bpic_recovery.csv")); print(R, row.names = FALSE)
