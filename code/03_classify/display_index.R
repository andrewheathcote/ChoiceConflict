## The paper figures' convergence classes (Andrew, 2026-10-20). A model is DISPLAYED AS CONVERGED when a current
## joint-move-family fit (jm, its Gadi continuations, the jm1 seeds, the AWS campaign jmh fits) either
##   - met the stop rule (max(mu Rhat, alpha .99 quantile) < 1.1), or
##   - is partially converged (called ridge-only before 2026-11-10: rhat_ident <= 1.1, LL drift <= 0.2) with the rule's Rhat < 1.5, not split within the run
##     (subject spread <= 15, chain LL gap <= 100) and without a subject-mode deficit (jm screen flag, or > 10 LL below the
##     best run of its case).
## Converged fits with a subject-mode caveat (split, or stuck in a subject's worse mode) are NOT displayed as converged
## (Andrew, 2026-10-20); neither are partially converged fits at Rhat >= 1.5.
## Writes derived/tables/display_jm_base.csv (03_classify/apply_q6.R then adds the Q6 refits): one row per model, its class, the run whose DIC the figures should use, and that DIC
## (jm1 seeds: from jm1_dic.csv; NA where a run's DIC has not been computed).
##   Rscript code/03_classify/display_index.R
## Inputs are the frozen outputs of the earlier steps (derived/tables/frozen/): conv_by_spec.csv (conv_by_sampler.R, dic_trust.R),
## inventory_classes.csv (replay/inventory_classes.R: graded runs, LL drift, ridge class = partially converged here), index_jm.csv (the jm index),
## screen_flags.csv and retry_eval.csv (the subject-mode screen), campaign_summary_dic.csv (the DIC of each campaign, Niek
## and local fit summary), jm1_dic.csv, handstopped_runs.csv, ll_screen_subjects.csv (RDM likelihood-defect screen).
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
fz <- function(f) read.csv(tab_path(file.path("frozen", f)), stringsAsFactors = FALSE)
U <- fz("conv_by_spec.csv"); A <- fz("inventory_classes.csv")
J <- fz("index_jm.csv"); F <- fz("screen_flags.csv")
## a summary's label is the run tag ("jmh1") or, in the full Niek summaries written by hardcases_niek/summarise_niek.R,
## the whole label ("ANT_RDM_BiVc_niek"): both give <spec>_<tag> here (2026-10-28)
H <- fz("campaign_summary_dic.csv")
## every current run with what is known about it
inv <- data.frame(spec = A$spec, label = A$label, source = A$source, met = A$met %in% TRUE, ridge = A$class == "ridge-only",
                  rule = pmax(A$rhat_mu, A$a_q99), split = A$subj_spread_max > 15 | A$gap_ll > 100, deficit = A$deficit)
k <- inv$source == "jm" & inv$spec %in% J$spec[J$label == "jm" & J$stop_rule %in% TRUE]; inv$met <- inv$met | k
jm <- data.frame(spec = J$spec, label = paste0(J$spec, "_", J$label), source = ifelse(J$label == "jm", "jm", "jm (continued)"),
                 met = J$stop_rule %in% TRUE, ridge = FALSE, rule = pmax(J$rhat_mu, J$rhat_alpha_q99), split = FALSE, deficit = NA)
jm <- jm[!(jm$source == "jm" & paste(jm$spec, "jm") %in% paste(inv$spec, inv$source)), ]
R <- rbind(inv, jm)
## the jm screen flags (screen_flags.csv) are about the original jm runs; the retries refit for the s348 / s454 trap
## (jm2-jm4) carry their own check (retry_eval.csv: good_mode), so a retry is flagged only when it stayed in the worse mode
RE <- fz("retry_eval.csv")
R$flag <- (R$source == "jm" & R$spec %in% F$spec) |
          (R$source == "jm (continued)" & paste0(R$spec, "_", sub(".*_", "", R$label)) %in% paste0(RE$spec, "_", RE$label)[!RE$good_mode])
## runs whose posterior sits where the pre-26-Oct RDM C++ likelihood is wrong (hardcases_niek/REPORT.md section 10,
## screen/ll_screen_subjects.csv: a subject's log-likelihood at its posterior mean > 5 above the fixed build's): their
## DIC is spuriously low, so they are not clean (2026-10-28)
LDF <- tab_path("frozen/ll_screen_subjects.csv")
LD <- if (file.exists(LDF)) read.csv(LDF, stringsAsFactors = FALSE) else data.frame(label = character(), spec = character(), subject = character(), diff = numeric())
LD <- LD[LD$diff > 5, ]; R$flag <- R$flag | R$label %in% LD$label
R$clean <- !R$split & !(!is.na(R$deficit) & R$deficit > 10) & !R$flag
## DIC of each run: the index's selected run (jm, or its continuation) from index_jm; an inventory jm run only when it is
## the index's selected run; campaign runs from their summaries; jm1 seeds not computed (NA)
R$DIC <- NA_real_
idx <- R$source %in% c("jm", "jm (continued)") & R$label == paste0(R$spec, "_", J$label[match(R$spec, J$spec)])
R$DIC[idx] <- J$DIC[match(R$spec[idx], J$spec)]
h <- R$source %in% c("jmh", "jmh (w3, record)", "niek", "local"); R$DIC[h] <- H$DIC[match(R$label[h], H$label)]
## jm1 seeds: DICs computed from their fits on S3 (jm1_dic.csv, 2026-10-20) for the runs the display uses
J1 <- fz("jm1_dic.csv"); j1 <- R$source == "jm1"; R$DIC[j1] <- J1$DIC[match(R$label[j1], J1$label)]
## class of each run, best first; a model takes its best run's class (ties: a run with a known DIC, then the lower rule Rhat)
lev <- c("converged", "converged, subject-mode caveat", "partially converged", "partially converged, Rhat >= 1.5", "not converged")
R$class <- ifelse(R$met & R$clean, lev[1], ifelse(R$met, lev[2], ifelse(R$ridge & R$clean & R$rule < 1.5, lev[3],
           ifelse(R$ridge & R$clean, lev[4], lev[5]))))
## runs stopped by hand with no end record (handstopped_runs.csv: ANT_DDM_BciVcSc jmh1c6, DIC computed from its fit):
## not converged, with that DIC (2026-11-02)
HS <- fz("handstopped_runs.csv")
R <- rbind(R, data.frame(spec = HS$spec, label = HS$label, source = HS$source, met = FALSE, ridge = FALSE, rule = HS$rule,
                         split = NA, deficit = NA, flag = FALSE, clean = FALSE, DIC = HS$DIC, class = lev[5])[, names(R)])
R <- R[order(R$spec, match(R$class, lev), is.na(R$DIC), R$rule), ]; B <- R[!duplicated(R$spec), ]
D <- U[, c("spec", "task", "arch", "in_play")]
D <- merge(D, B[, c("spec", "class", "label", "source", "rule", "DIC")], by = "spec", all.x = TRUE)
D$class[is.na(D$class)] <- "not fitted yet"
D$class[D$in_play & !D$class %in% lev[c(1, 3)]] <- "in play (campaign fit queued or running)"
D$display_converged <- D$class %in% lev[c(1, 3)]                     # Andrew 2026-10-20: the subject-mode caveat models are not shown as converged
D$caveat <- ifelse(D$class == lev[2], paste("subject mode:", vapply(D$spec, function(s) paste(unique(F$subject[F$spec == s]), collapse = " "), "")), "")
D$caveat[D$class == lev[2] & D$caveat == "subject mode: "] <- "subject mode: split run or > 10 LL below another run of the case"
ld <- D$label %in% LD$label
D$caveat[ld] <- paste("likelihood defect (pre-26-Oct RDM C++ likelihood):", vapply(D$spec[ld], function(s) paste(unique(LD$subject[LD$spec == s]), collapse = " "), ""))
D <- D[order(D$task, D$arch, D$spec), ]
out <- tab_path("display_jm_base.csv"); write_out(D, out)
cat("wrote", normalizePath(out), "\n"); options(width = 200)
print(table(paste(D$task, D$arch), factor(D$class, c(lev, "in play (campaign fit queued or running)", "not fitted yet"), labels = c("conv", "conv*", "partial<1.5", "partial>=1.5", "not", "in play", "not fitted"))))
cat("\ndisplayed as converged:", sum(D$display_converged), "of", nrow(D), "; DIC missing for displayed models:", sum(D$display_converged & is.na(D$DIC)), ":",
    paste(D$spec[D$display_converged & is.na(D$DIC)], "(", D$label[D$display_converged & is.na(D$DIC)], ")", collapse = ", "), "\n")
