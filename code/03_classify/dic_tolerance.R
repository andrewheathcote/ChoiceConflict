## How far a partially converged fit's DIC can be from the converged answer (supplement S6). Both checks use the per-subject
## log-likelihoods every fit stores (samples$subj_ll), so no likelihood is re-evaluated and no EMC2 build is loaded:
## best-draw DIC = sum_s (2 meanD_s - minD_s) (pD_s = meanD_s - minD_s; the minimum-deviance draw is never worse than the
## posterior-mean point in these fits), which reproduces the display DIC exactly.
##   1. Pairs: every superseded run that is partially (or not) converged and whose model's display run met the stop rule
##      (FITS.csv rows with in_display FALSE). Both runs are thinned evenly within chains to the same total draw count;
##      the DIC difference (superseded minus converged; > 0: the unconverged run makes the model look worse) is split
##      into mean deviance and pD.
##   2. Chains: for every conflict fit in FITS.csv, DIC from each chain on its own, each chain cut evenly to the shortest
##      chain's draw count. The spread across chains shows how much DIC depends on where the chains sit; converged fits
##      give the Monte Carlo baseline. Diffusion models have no converged twin, so this is their evidence.
##   Rscript code/03_classify/dic_tolerance.R [cores]
## Reads derived/tables/{FITS,display_jm}.csv, frozen/inventory_classes.csv (rule value of superseded runs) and the fits;
## writes derived/tables/dic_tolerance_{pairs,chains}.csv. Per-fit results are cached in derived/dic_tolerance/.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
cores <- suppressWarnings(as.integer(commandArgs(TRUE)[1])); if (is.na(cores)) cores <- 4
CACHE <- file.path(DERIVED, "dic_tolerance"); dir.create(CACHE, showWarnings = FALSE)
F <- read_tab("FITS.csv"); F <- F[F$task %in% TASKS, ]
D <- read_tab("display_jm.csv")
I <- read.csv(tab_path(file.path("frozen", "inventory_classes.csv")), stringsAsFactors = FALSE)
F$rule <- ifelse(F$in_display, D$rule[match(F$label, D$label)], pmax(I$rhat_mu, I$a_q99)[match(F$label, I$label)])

## sample-stage log-likelihoods, one subjects x draws matrix per chain (the fit object is a plain list; no package needed)
fit_ll <- function(label) { f <- load_fit(label); lapply(f, function(ch) ch$samples$subj_ll[, ch$samples$stage == "sample", drop = FALSE]) }
dic_of <- function(L) { Dv <- -2 * do.call(cbind, L); mD <- rowMeans(Dv); mn <- apply(Dv, 1, min)
  c(DIC = sum(2 * mD - mn), meanD = sum(mD), pD = sum(mD - mn)) }
## evenly spaced draws within each chain, as many in total as `total`, shared out in proportion to the chains' lengths
thin <- function(L, total) { n <- sapply(L, ncol); k <- floor(total * n / sum(n)); r <- total - sum(k)
  if (r > 0) { o <- order(n - k, decreasing = TRUE)[seq_len(r)]; k[o] <- k[o] + 1 }
  mapply(function(x, m) x[, round(seq(1, ncol(x), length.out = m)), drop = FALSE], L, k, SIMPLIFY = FALSE) }
chain_stats <- function(label, L) {
  f <- file.path(CACHE, paste0(label, ".rds")); if (file.exists(f)) return(readRDS(f))
  m <- min(sapply(L, ncol)); C <- sapply(L, function(x) dic_of(list(thin(list(x), m)[[1]])))
  s <- list(label = label, chains = length(L), draws = sum(sapply(L, ncol)), draws_chain = m, full = dic_of(L), per_chain = C)
  saveRDS(s, f); s }

## 1. pairs
S <- F[!F$in_display, ]
S$partner <- D$label[match(S$spec, D$spec)]; S$partner_class <- D$class[match(S$spec, D$spec)]
S <- S[S$partner_class == "converged", ]
pair <- function(i) {
  f <- file.path(CACHE, paste0("pair_", S$label[i], ".rds")); if (file.exists(f)) return(readRDS(f))
  L1 <- fit_ll(S$label[i]); a <- chain_stats(S$label[i], L1)
  L2 <- fit_ll(S$partner[i]); b <- chain_stats(S$partner[i], L2)
  tot <- min(a$draws, b$draws); x <- dic_of(thin(L1, tot)); y <- dic_of(thin(L2, tot))
  r <- data.frame(spec = S$spec[i], task = S$task[i], arch = S$arch[i], run = S$label[i], class = S$class[i],
                  rule = round(S$rule[i], 3), converged_run = S$partner[i], draws = tot,
                  DIC_run = x[["DIC"]], DIC_converged = y[["DIC"]], diff_DIC = x[["DIC"]] - y[["DIC"]],
                  diff_meanD = x[["meanD"]] - y[["meanD"]], diff_pD = x[["pD"]] - y[["pD"]],
                  run_chain_sd = sd(a$per_chain["DIC", ]), converged_chain_sd = sd(b$per_chain["DIC", ]))
  saveRDS(r, f); r }
P <- do.call(rbind, parallel::mclapply(seq_len(nrow(S)), pair, mc.cores = cores, mc.preschedule = FALSE))
write_out(P, tab_path("dic_tolerance_pairs.csv"))

## 2. chains, for every conflict fit (the pairs' fits come from the cache)
one <- function(label) { s <- if (file.exists(file.path(CACHE, paste0(label, ".rds")))) readRDS(file.path(CACHE, paste0(label, ".rds")))
                         else chain_stats(label, fit_ll(label))
  C <- s$per_chain
  data.frame(label = label, chains = s$chains, draws = s$draws, draws_chain = s$draws_chain, DIC = s$full[["DIC"]],
             chain_sd_DIC = sd(C["DIC", ]), chain_range_DIC = diff(range(C["DIC", ])),
             chain_sd_meanD = sd(C["meanD", ]), chain_sd_pD = sd(C["pD", ])) }
big <- order(-F$size)                                                    # largest first, so the tail is short
G <- do.call(rbind, parallel::mclapply(F$label[big], one, mc.cores = cores, mc.preschedule = FALSE))
G <- merge(F[, c("label", "spec", "task", "arch", "class", "in_display", "rule")], G, by = "label")
G$display_DIC <- D$DIC[match(G$label, D$label)]                          # display DIC reproduced?
cat("display DIC reproduced (within 0.01):", sum(abs(G$DIC - G$display_DIC) < 0.01, na.rm = TRUE), "of", sum(!is.na(G$display_DIC)), "\n")
x <- G[!is.na(G$display_DIC) & abs(G$DIC - G$display_DIC) >= 0.01, ]; if (nrow(x)) print(x[, c("label", "DIC", "display_DIC")], row.names = FALSE)
write_out(G[order(G$task, G$arch, G$label), ], tab_path("dic_tolerance_chains.csv"))

## summaries
fam <- function(a) ifelse(a %in% c("DDM", "WDM"), "diffusion", "race")
cat("\n=== 1. Pairs: run minus converged run of the same model, at matched draws\n")
P$band <- cut(P$rule, c(1, 1.3, 1.5, Inf), right = FALSE, labels = c("1.1-1.3", "1.3-1.5", ">= 1.5"))
for (g in list(P$arch, P$band, P$class)) print(do.call(rbind, lapply(split(P, g, drop = TRUE), function(x)
  data.frame(n = nrow(x), median = round(median(x$diff_DIC)), max_abs = round(max(abs(x$diff_DIC))),
             max_abs_meanD = round(max(abs(x$diff_meanD))), max_abs_pD = round(max(abs(x$diff_pD)))))))
cat("\n|diff_DIC| / the run's chain SD: max", round(max(abs(P$diff_DIC) / P$run_chain_sd), 2),
    " median", round(median(abs(P$diff_DIC) / P$run_chain_sd), 2), "\n")
cat("\n=== 2. Chain SD of DIC, displayed fits, by class and family\n")
H <- G[G$in_display, ]
print(do.call(rbind, lapply(split(H, list(H$class, fam(H$arch)), drop = TRUE), function(x)
  data.frame(n = nrow(x), median_sd = round(median(x$chain_sd_DIC)), max_sd = round(max(x$chain_sd_DIC)),
             median_sd_meanD = round(median(x$chain_sd_meanD)), median_sd_pD = round(median(x$chain_sd_pD))))))
