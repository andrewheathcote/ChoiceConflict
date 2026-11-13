## PSIS-LOO (and WAIC from the same pointwise values) for every real conflict fit in derived/tables/display_jm.csv,
## at n_draws = 2000, on DMC2 f2f58985 (lib/f2f58985). Each fit runs as its own process
## (code/08_recovery/conflict/loo.R); this script schedules them under a core and memory budget and collates.
##   Rscript code/04_select/loo_real.R [cores = 8] [mem_GB = 60] [ddm_threads = 2]
##   Rscript code/04_select/loo_real.R collate
## Inputs: derived/tables/display_jm.csv, fits/<task>/<arch>/<label>.RData (read only).
## Outputs: derived/loo_real/<label>.rds (per fit, with pointwise elpd), derived/loo_real/logs/<label>.log,
##          derived/tables/loo_real.csv (collate; one row per fit with spec, task, arch and DIC from display_jm.csv).
## Resumable: a fit whose .rds exists, or that is listed in derived/loo_real/elsewhere.txt, is skipped. Memory per fit is budgeted as 2 + 12 x its .RData size (GB), which
## matched the measured peaks (work/loo_full/TIMING.md). DDM (diffusion with the slow likelihood) runs first.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
a <- commandArgs(TRUE)
OUT <- file.path(DERIVED, "loo_real"); LOGS <- file.path(OUT, "logs"); dir.create(LOGS, recursive = TRUE, showWarnings = FALSE)
D <- read_tab("display_jm.csv")
D$file <- vapply(D$label, fit_path, "")
collate <- function() {
  f <- file.path(OUT, paste0(D$label, ".rds")); ok <- file.exists(f)
  S <- do.call(rbind, lapply(f[ok], function(x) readRDS(x)$summary))
  R <- merge(D[ok, c("spec", "task", "arch", "label", "DIC", "display_converged")], S, by = "label", sort = FALSE)
  write_out(R[order(R$task, R$arch, R$LOO), ], tab_path("loo_real.csv"))
  cat(sum(ok), "of", nrow(D), "fits collated into", tab_path("loo_real.csv"), "\n")
}
if (length(a) && a[1] == "collate") { collate(); quit(save = "no") }
CORES <- if (length(a) >= 1) as.integer(a[1]) else 8L
MEM   <- if (length(a) >= 2) as.numeric(a[2]) else 60
DDMT  <- if (length(a) >= 3) as.integer(a[3]) else 2L
LIB <- file.path(LIBDIR, "f2f58985"); stopifnot(startsWith(readLines(file.path(LIB, "EMC2_SHA.txt"))[1], "f2f58985"))
RUN <- file.path(CODE, "08_recovery", "conflict", "loo.R")
TIME <- if (Sys.info()[["sysname"]] == "Darwin") "/usr/bin/time -l" else "/usr/bin/time -v"
D$mem <- 2 + ifelse(D$arch == "DDM", 9, 12) * file.size(D$file) / 1e9   # DDM: 8.6 GB measured for a 0.73 GB ANT fit
D$threads <- ifelse(D$arch == "DDM", DDMT, 1L)
## fits run elsewhere (derived/loo_real/elsewhere.txt, one label per line; their .rds are copied in afterwards), and
## loo.R processes already running (e.g. from an earlier driver), whose cores and memory count against the budget
skip <- if (file.exists(file.path(OUT, "elsewhere.txt"))) readLines(file.path(OUT, "elsewhere.txt")) else character(0)
ext_running <- function() { p <- system("ps -axo command", intern = TRUE); p <- p[grepl("conflict/loo.R", p) & grepl("exec/R", p)]
  which(vapply(D$label, function(l) any(grepl(paste0("/", l, ".rds"), p, fixed = TRUE)), TRUE)) }
## DDM largest first, then the rest alternating largest / smallest so the memory budget packs
o1 <- which(D$arch == "DDM"); o1 <- o1[order(-D$mem[o1])]
o2 <- which(D$arch != "DDM"); o2 <- o2[order(-D$mem[o2])]
alt <- as.vector(rbind(o2, rev(o2)))[seq_along(o2)]; alt <- unique(c(alt, o2))
Q <- c(o1, alt)
Q <- Q[!file.exists(file.path(OUT, paste0(D$label[Q], ".rds"))) & !D$label[Q] %in% skip]
Q <- setdiff(Q, ext_running())
cat(sprintf("%d fits to run (%d done), %d cores, %.0f GB, DDM threads %d\n", length(Q), nrow(D) - length(Q), CORES, MEM, DDMT))
running <- integer(0); failed <- character(0); t0 <- Sys.time()
exitf <- function(i) file.path(LOGS, paste0(D$label[i], ".exit"))
while (length(Q) || length(running)) {
  for (i in running) if (file.exists(exitf(i))) {
    st <- readLines(exitf(i)); running <- setdiff(running, i)
    if (!identical(st, "0") || !file.exists(file.path(OUT, paste0(D$label[i], ".rds")))) failed <- c(failed, D$label[i])
    cat(format(Sys.time(), "%H:%M"), D$label[i], if (identical(st, "0")) "done" else paste("FAILED", st),
        sprintf("(%d left, %d running)\n", length(Q), length(running)))
  }
  ext <- setdiff(ext_running(), running)
  repeat {
    if (!length(Q)) break
    busy <- union(running, ext); i <- Q[1]; uc <- sum(D$threads[busy]); um <- sum(D$mem[busy])
    if (length(busy) &&(uc + D$threads[i] > CORES || um + D$mem[i] > MEM)) break
    Q <- Q[-1]; unlink(exitf(i)); running <- c(running, i)
    cmd <- sprintf("EMC2_LIB=%s OMP_NUM_THREADS=%d %s Rscript %s %s %s 2000 %d > %s 2>&1; echo $? > %s",
                   shQuote(LIB), D$threads[i], TIME, shQuote(RUN), shQuote(D$file[i]),
                   shQuote(file.path(OUT, paste0(D$label[i], ".rds"))), D$threads[i],
                   shQuote(file.path(LOGS, paste0(D$label[i], ".log"))), shQuote(exitf(i)))
    system(paste("sh -c", shQuote(cmd)), wait = FALSE)
  }
  Sys.sleep(5)
}
cat(sprintf("all done in %.0f min; %d failed%s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), length(failed),
            if (length(failed)) paste0(": ", paste(failed, collapse = " ")) else ""))
collate()
