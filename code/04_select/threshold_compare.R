## Directional DIC thresholds against the single floor, from derived/tables/display_jm.csv. Read-only on the fits.
##   Rscript code/04_select/threshold_compare.R     -> derived/tables/threshold_compare_{pairs,mechs}.csv
## Single floor (05_floors/recovery_floor.R): any difference > 150 counts (the largest recovery gap, 149.4, rounded up). Directional (recovery studies 1-5, 3 Nov): a model with MORE
## mechanisms is supported over one with fewer if it leads by > ADD (50; complex wrong winners reached 40 and 28); a
## model with FEWER mechanisms, or a same-size rival, or another architecture, needs > DROP (149; simpler wrong winners
## reached 149.4 and 53; DROP is the floor, 150). Mechanisms are the tokens of floor_within_jm.R (Vci = Vc + Vi). Only displayed models
## (display_converged); the best of an architecture may be partially converged (Andrew, 2026-11-10; it had to meet the
## stop rule before), as in the heat grids; (p) marks partially converged.
ADD <- 50; DROP <- 150
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
d <- read_tab("display_jm.csv"); d <- d[!is.na(d$DIC) & d$display_converged, ]
d$var <- sub("^[^_]+_[^_]+_", "", d$spec)
tok <- function(v) { m <- regmatches(v, gregexpr("(SV|[A-Z])[ci]*", v))[[1]]
  unlist(lapply(m, function(t) { p <- sub("[ci]+$", "", t); s <- sub("^(SV|[A-Z])", "", t); paste0(p, strsplit(s, "")[[1]]) })) }
d$tk <- lapply(d$var, tok); d$n <- lengths(d$tk); d$met <- d$class == "converged"
lab <- function(x) paste0(x$var, ifelse(x$met, "", "(p)"))
TA <- expand.grid(arch = c("LNR", "RDM", "LBA", "WDM", "DDM"), task = c("Simon", "Stroop", "ANT"), stringsAsFactors = FALSE)[, 2:1]

cat("=== 1. Nested pairs: child = parent + one mechanism. Changes are pairs where the child leads by", ADD, "to", DROP, "\n")
P <- NULL
for (r in seq_len(nrow(TA))) { x <- d[d$task == TA$task[r] & d$arch == TA$arch[r], ]
  for (i in seq_len(nrow(x))) for (j in seq_len(nrow(x))) {
    if (x$n[j] != x$n[i] + 1 || !all(x$tk[[i]] %in% x$tk[[j]])) next
    g <- x$DIC[i] - x$DIC[j]                                                   # > 0: the child (more mechanisms) leads
    P <- rbind(P, data.frame(task = x$task[i], arch = x$arch[i], parent = lab(x[i, ]), child = lab(x[j, ]),
                             added = setdiff(x$tk[[j]], x$tk[[i]]), child_lead = round(g),
                             single = if (g > DROP) "added supported" else if (-g > DROP) "not needed" else "unresolved",
                             directional = if (g > ADD) "added supported" else if (-g > DROP) "not needed" else "unresolved")) } }
cat("nested pairs:", nrow(P), "\n"); print(table(single = P$single, directional = P$directional))
ch <- P[P$single != P$directional, ]; ch <- ch[order(ch$task, ch$arch, -ch$child_lead), ]
cat("\nnewly supported additions (child leads by", ADD, "-", DROP, "):\n"); print(ch[, c("task", "arch", "parent", "child", "added", "child_lead")], row.names = FALSE)
write_out(P, tab_path("threshold_compare_pairs.csv"))

cat("\n=== 2. Mechanism necessity: best model without m vs the architecture's best (met rule).\n")
cat("    single: needed if the gap > 150. directional: > 50 when the best model without m has fewer mechanisms than the best, else > 150\n")
N <- NULL
for (r in seq_len(nrow(TA))) { x <- d[d$task == TA$task[r] & d$arch == TA$arch[r], ]; if (!nrow(x)) next
  b <- x[which.min(x$DIC), ]
  for (m in b$tk[[1]]) { wo <- x[!sapply(x$tk, function(t) m %in% t), ]; if (!nrow(wo)) next
    w <- wo[which.min(wo$DIC), ]; g <- w$DIC - b$DIC; thr <- if (w$n < b$n) ADD else DROP
    N <- rbind(N, data.frame(task = b$task, arch = b$arch, best = lab(b), mech = m, best_without = lab(w), gap = round(g),
                             fewer = w$n < b$n, single = g > DROP, directional = g > thr)) } }
print(N[N$single | N$directional, ], row.names = FALSE)
cat("\nchanged by the directional rule:\n"); print(N[N$single != N$directional, ], row.names = FALSE)
write_out(N, tab_path("threshold_compare_mechs.csv"))

cat("\n=== 3. Contenders within an architecture: models not separable from its best (met rule).\n")
cat("    single: within 150. directional: within 50 if it has fewer mechanisms than the best, else within 150\n")
C <- NULL
for (r in seq_len(nrow(TA))) { x <- d[d$task == TA$task[r] & d$arch == TA$arch[r], ]; if (!nrow(x)) next
  b <- x[which.min(x$DIC), ]; o <- x[x$spec != b$spec, ]; g <- o$DIC - b$DIC
  fewer <- o$n < b$n                                                        # fewer mechanisms, nested or not (study 5: Bi vs ViSi)
  s1 <- g <= DROP; s2 <- g <= ifelse(fewer, ADD, DROP)
  C <- rbind(C, data.frame(task = b$task, arch = b$arch, best = lab(b), shown = nrow(x), single = sum(s1), directional = sum(s2),
                           dropped = paste(sprintf("%s+%.0f", lab(o[s1 & !s2, ]), g[s1 & !s2]), collapse = " "))) }
print(C, row.names = FALSE)
