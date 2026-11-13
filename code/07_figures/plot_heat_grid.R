## Combined figure: a 4-column x 2-row grid.
##   top    LNR | WDM | DDM | legend        (four-mechanism, 1-4-6-4 each)
##   bottom RDM (2 cols) | LBA (2 cols)     (three two-parameter trees each)
## Every panel shares one label size, computed from the tightest level anywhere.
## Values are uniform random on 0-100, standing in for a real misfit metric.
##
## DISPLAY VERSION (2026-10-20): copy of plot_heat_grid_jm.R.
##   [PP_AGG=median] Rscript code/07_figures/plot_heat_grid.R dic|delta|caf
## PP_AGG (delta/caf only): group curve = mean (default) or median over subjects; see pp_display_lib.R.
## Convergence and DIC come from derived/tables/display_jm.csv (03_classify). Delta and CAF misfit come from
## derived/tables/adequacy_display[_median].csv (06_predict/build_adequacy.R: each displayed fit's SUBJECT-LEVEL
## posterior predictive).
## Cell states: converged; partially converged (stop rule not met, Rhat_ident <= 1.1, rule Rhat < 1.5) (same fill, hatched);
## did not converge (grey: not converged, partially converged with Rhat >= 1.5, and the
## subject-mode-caveat models); fit running (campaign DDMs in play, dotted).
## mech comes from derived/tables/mech_key.csv by spec.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
source(file.path(CODE, "lib", "name_scheme.R"))
MEASURE <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "dic"
## Joint DIC is dropped from this project: it recovered the generating model
## 0 of 20 times in the RDM recovery study and 2 of 20 (below the 20% chance
## rate) in the LNR one, so it is not a weak criterion but a misleading one.
## ic_joint.csv and the code that built it are kept as the evidence for that
## decision; nothing downstream reports it.
stopifnot(MEASURE %in% c("delta","caf","dic"))
if (MEASURE == "jdic") stop("joint DIC has been dropped from this project")
AGG <- Sys.getenv("PP_AGG", unset = "mean"); stopifnot(AGG %in% c("mean", "median"))
aggsuf <- if (AGG == "mean" || MEASURE == "dic") "" else "_median"
## May a partially converged model take the best box, the best-2 circle, the zero point and
## the model-type header value? TRUE (Andrew, 2026-11-10; FALSE from 2026-10-20): their DIC is good to
## about +-75 (03_classify/dic_tolerance.R), and on ANT the two best DDMs are partially converged, 1,221 and 582 better than the best
## DDM that met the stop rule. They keep their hatching.
RIDGE_CAN_BE_BEST <- TRUE
TASKS <- strsplit(Sys.getenv("HEAT_TASKS", unset = "Simon,Stroop,ANT"), ",")[[1]]   # HEAT_TASKS=Simon for a one-page draft
OUTDIR <- Sys.getenv("HEAT_OUT", unset = FIGDIR)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(OUTDIR, sprintf("heat_grid_%s.pdf", MEASURE))
FOUR  <- list(WDM = c("a","z","v","s"), DDM = c("a","z","v","s"), LNR = c("M","m","S","s"))
ROOTS <- list(VS = c("V","v","S","s"), BS = c("B","b","S","s"), BV = c("B","b","V","v"))
KEEP  <- c(VS = "Vci", BS = "Sci", BV = "Bci")
BPAR  <- c(V = "VcSi", v = "ViSc", S = "Sci", s = "BcSi", B = "Bci", b = "BiVi")
SWAP  <- list(c("VcSi", "ViSi"))
BH <- 0.020; MAXCEX <- 1.05
## Reversed: index 1 is the BEST value, so the best cells are darkest red and
## the worst darkest blue. hcl.colors gives NA at an odd-n midpoint, hence the
## explicit ramp.
RAMP <- colorRampPalette(rev(c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8",
                           "#FFFFBF","#FEE090","#FDAE61","#F46D43","#D73027",
                           "#A50026")))(101)
stopifnot(!any(is.na(RAMP)))
## ---- values, one table per task -------------------------------------------
## delta / caf : adequacy statistics for that task
## dic / jdic  : conditional or joint DIC as a difference from the best model in
##               that task (0 = best). The colour scale is shared across the
##               three pages so tasks are comparable, and capped at the 98th
##               percentile so one runaway cell cannot flatten the rest.
skey <- function(m) sapply(strsplit(m, ""), function(z) paste(sort(z), collapse = ""))
DISP <- read_tab("display_jm.csv")
local({
  MK <- read_tab("mech_key.csv")              # one mechanism key (03_classify/mech_key.R)
  DISP$variant <- sub("^[^_]+_", "", DISP$spec)
  DISP$mech <- MK$mech[match(DISP$spec, MK$spec)]
  stopifnot(!anyNA(DISP$mech), all(nzchar(DISP$mech)))
  DISP$ridge <- DISP$class == "partially converged"
  DISP$kk <- paste(DISP$task, DISP$arch, skey(DISP$mech))
  stopifnot(!anyDuplicated(DISP$kk))
  AD <- read_tab(sprintf("adequacy_display%s.csv", aggsuf))
  stopifnot(!anyDuplicated(AD$spec))
  for (v in c("rel","caf","degen")) DISP[[v]] <- setNames(AD[[v]], AD$spec)[DISP$spec]
  DISP <<- DISP })
vals_for <- function(tk) {
  if (MEASURE %in% c("delta","caf")) {
    ## displayed-converged models only, partially converged included (they are never marked best when
    ## RIDGE_CAN_BE_BEST is FALSE: nobest() below); misfit has no zero point to protect
    D <- DISP
    D$v <- if (MEASURE == "delta") D$rel else D$caf
    D <- D[!is.na(D$v) & D$task == tk & D$display_converged & D$degen == 0, ]
    D <- D[!(D$arch %in% c("WDM","DDM") & grepl("V", D$mech)), ]
    setNames(D$v, paste(D$arch, skey(D$mech)))
  } else {
    ## joint-move family, best current fit per model: figs/display_jm.csv
    D <- DISP
    D$v <- D$DIC
    ## display_converged replaces the stop_rule test (same role: a model not displayed as
    ## converged must not set the zero point)
    D <- D[!is.na(D$v) & D$task == tk & D$display_converged, ]
    if (!RIDGE_CAN_BE_BEST) D <- D[!D$ridge, ]
    D <- D[!(D$arch %in% c("WDM","DDM") & grepl("V", D$mech)), ]
    D$val <- D$v - min(D$v); D$kk <- paste(D$arch, skey(D$mech))
    A <- aggregate(val ~ kk, D, min); setNames(A$val, A$kk)
  }
}
## partially converged DIC (good to about +-75: 03_classify/dic_tolerance.R) is added after the reference is set, so
## with RIDGE_CAN_BE_BEST = FALSE it is coloured but never the zero point
ridge_vals <- function(tk, ref) {
  D <- DISP[DISP$ridge & DISP$task == tk & !is.na(DISP$DIC), ]
  D <- D[!(D$arch %in% c("WDM","DDM") & grepl("V", D$mech)), ]
  setNames(D$DIC - ref, paste(D$arch, skey(D$mech)))
}
VL     <- lapply(TASKS, vals_for); names(VL) <- TASKS
if (!RIDGE_CAN_BE_BEST && MEASURE == "dic") for (tk in TASKS) {
  ref <- min(DISP$DIC[DISP$task == tk & DISP$display_converged & !DISP$ridge], na.rm = TRUE)
  VL[[tk]] <- c(VL[[tk]], ridge_vals(tk, ref)) }
NFL <- read.csv(file.path(INPUTS, "noise_floors.csv"), stringsAsFactors = FALSE)   # DIC floor 150: the largest recovery gap, 149.4, rounded up (05_floors/recovery_floor.R)
## delta / CAF floors for the page's own aggregation (mean or median group curve; noise_floor_agg.R, Andrew 4 Nov):
## the median delta floor is about twice the mean one. DIC keeps noise_floors.csv.
NFA <- read_tab("noise_floors_agg.csv")
nfloor <- function(meas, tk) {
  if (meas %in% c("delta", "caf")) { v <- NFA$floor[NFA$measure == meas & NFA$task == tk & NFA$agg == AGG]; return(if (length(v)) v[1] else NA_real_) }
  if (is.null(NFL)) return(NA_real_)
  v <- NFL$floor[NFL$measure == meas & NFL$task == tk]; if (length(v)) v[1] else NA_real_ }

## LOGOFF is reset per page below. For DIC the log offset is the difference
## floor, not 1: with an offset of 1 the meaningless region below the floor
## swallowed most of the scale, leaving ANT 97% in one half of the ramp.
LOGOFF <- if (MEASURE %in% c("delta","caf")) 0 else 1
allv   <- unlist(VL, use.names = FALSE)
## Global limits are a fallback only: LO/HI are recomputed per page below,
## so each task spans its own min to max and every page uses the whole ramp.
## A pooled scale let Stroop dominate -- its CAF misfit reaches 30.8 against
## Simon's 3.9, so the Simon and ANT CAF pages never got past 63 of 100.
LO <- min(allv); HI <- max(allv)
tf <- function(v) log10(pmax(v, LO) + LOGOFF)
## Outcomes per cell:
##   a value      displayed as converged (hatched: partially converged)
##   grey fill    fitted, not displayed as converged
##   dashed box   campaign fit still running
##   dotted box   no joint-move fit (Simon DDM_BciVcSc, never refit)
##   crossed box  NOT APPLICABLE by design -- only Stroop WDM/DDM variants that
##                carry the start point, which accuracy coding removes
## ANT DDM_BiVcSc is the one applicable cell never fitted; it is being fitted
## now and is shown as "did not converge" in the meantime.
FITTED  <- DISP$kk[!DISP$in_play]
RUNNING <- DISP$kk[DISP$in_play]
RIDGE   <- DISP$kk[DISP$ridge]
## The set of cells the figure actually draws. BEST must be chosen from these:
## the data also hold models outside the design (k > 4 race models, now in
## largek/), and on conditional DIC for Simon the overall minimum is one of
## those, so an unrestricted argmin marked a cell that is never drawn.
.lat <- function(L, n) unlist(lapply(1:n, function(k)
          apply(combn(L, k), 2, function(z) paste(sort(z), collapse = ""))))
.PR <- list(B=c("B","b"), V=c("V","v"), S=c("S","s"))
.np <- function(m) { l <- strsplit(m,"")[[1]]; sum(sapply(.PR, function(p) sum(p %in% l) == 1)) }
.r6 <- .lat(c("B","b","V","v","S","s"), 4)
.r6 <- .r6[sapply(.r6, function(m) nchar(m) <= 2 || .np(m) <= 1)]
DRAWN <- c(paste("WDM", .lat(c("a","z","v","s"), 4)),
           paste("DDM", .lat(c("a","z","v","s"), 4)),
           paste("LNR", .lat(c("M","m","S","s"), 4)),
           paste("RDM", .r6), paste("LBA", .r6))
VAL <- NULL; TK <- NULL; BEST <- NULL; BEST2 <- NULL   # set per page
NFCELL <- logical(0)
## TRUE for value keys ("arch mech") that may not be marked best on this page
nobest <- function(keys) if (RIDGE_CAN_BE_BEST) rep(FALSE, length(keys)) else paste(TK, keys) %in% RIDGE
getv <- function(arch, set) unname(VAL[paste(arch, paste(sort(set), collapse = ""))])
isbest <- function(arch, set)
  identical(paste(arch, paste(sort(set), collapse = "")), BEST)
getst <- function(arch, set) {
  k <- paste(TK, arch, paste(sort(set), collapse = ""))
  if (!is.na(getv(arch, set))) return(if (k %in% RIDGE) "rg" else "ok")
  if (k %in% RUNNING) return("rn")
  if (k %in% FITTED) return("nc")
  ## no jm fit exists for this cell: drawn as did not converge and reported below (Andrew, 3 Nov; was its own dotted
  ## "no jm fit" state). None since Simon DDM_BciVcSc was refit on Gadi (4 Nov, jm, not converged: it is in FITTED)
  if (!(TK == "Stroop" && arch %in% c("WDM","DDM") && "z" %in% set)) { NFCELL[k] <<- TRUE; return("nc") }
  ## applicable-by-design test, stated rather than inferred from absence
  if (TK == "Stroop" && arch %in% c("WDM","DDM") && "z" %in% set) "na" else "nc"
}
## Each model type is labelled with the value of its own best drawn cell --
## for RDM and LBA that is the minimum over all three of their trees, matching
## where the thin box goes.
archmin <- function(arch) {
  vd <- VAL[names(VAL) %in% DRAWN]
  vd <- vd[!nobest(names(vd))]
  v  <- vd[startsWith(names(vd), paste0(arch, " "))]
  if (!length(v) || all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
}
fmtv <- function(v) if (is.na(v)) "-" else switch(MEASURE,
  dic   = formatC(round(v), format = "d", big.mark = ","),   # "d" truncates: 75.9 -> 75
  jdic  = formatC(round(v), format = "d", big.mark = ","),
  delta = formatC(v, format = "f", digits = 1),
          formatC(v, format = "f", digits = 2))
## "(=)" after a model type whose best is within 150 of the page's best (not the best type itself), as the cell rule
## delta / CAF (4 Nov): within the page's difference floor PNF of the page's best misfit LO
archlab <- function(arch) { v <- archmin(arch)
  sep <- if (is.na(v)) FALSE else if (MEASURE == "dic") v > 0 && v <= DROP_THR else v > LO && v - LO <= PNF
  sprintf("%s: %s%s", arch, fmtv(v), if (isTRUE(sep)) " (=)" else "") }
col_of <- function(v) if (is.na(v)) "grey88" else
  RAMP[pmax(1, pmin(101, round(100*(tf(pmin(v,HI)) - tf(LO))/(tf(HI) - tf(LO))) + 1))]
cell <- function(xc, y, w, v, st = "ok", best = FALSE) {
  x1 <- xc-w/2*0.94; x2 <- xc+w/2*0.94; y1 <- y-BH/2; y2 <- y+BH/2
  if (st == "na") {                            # not applicable: crossed, unfilled
    rect(x1, y1, x2, y2, col = NA, border = "grey55", lwd = 0.5)
    segments(x1, y1, x2, y2, col = "grey55", lwd = 0.5)
    segments(x1, y2, x2, y1, col = "grey55", lwd = 0.5)
  } else if (st == "nf") {                     # no jm fit
    rect(x1, y1, x2, y2, col = NA, border = "grey70", lwd = 0.5, lty = 3)
  } else if (st == "rn") {                     # campaign fit queued or running
    rect(x1, y1, x2, y2, col = "grey95", border = "grey50", lwd = 0.5, lty = 2)
  } else if (st == "nc") {                     # fitted, did not converge
    rect(x1, y1, x2, y2, col = "grey82", border = "grey62", lwd = 0.4)
  } else { rect(x1, y1, x2, y2, col = col_of(v), border = NA)
    if (st == "rg") ridgehatch(x1, y1, x2, y2) }   # partially converged
  STATES <<- union(STATES, st)
}
## partially converged: slanted lines over the colour (was a white dot; the cell centre now carries the "=" mark, 3 Nov)
ridgehatch <- function(x1, y1, x2, y2) rect(x1, y1, x2, y2, density = 22, angle = 45, col = "grey15", lwd = 0.5, border = NA)
STATES <- character(0)
## the best cell is boxed together with its label, which makes it far easier to
## pick out than a border on the bar alone
## thick = best anywhere on the page, thin = best within its own tree
bestbox <- function(xc, ycell, ylab, w, lh, lwd = 3.8) {
  ## the page-best box (thick) is drawn larger than the within-tree boxes, and
  ## padded so its border clears the label rather than running against it
  big <- lwd > 2
  hw  <- max(w/2*0.94, 0.30*w) * (if (big) 1.22 else 1)
  pad <- if (big) lh*0.42 else 0
  yl  <- min(ycell - BH/2, ylab - lh*0.75) - pad
  yh  <- max(ycell + BH/2, ylab + lh*0.75) + pad
  rect(xc - hw, yl, xc + hw, yh, col = NA, border = "black", lwd = lwd)
}
## Best two-mechanism model, marked with an ellipse rather than a box so the
## two annotations never read as the same thing. Same size convention as
## bestbox: thick and larger for the best across all model types on the page,
## thin and tighter for the best within one model type.
bestcircle <- function(xc, ycell, ylab, w, lh, lwd = 3.8) {
  big <- lwd > 2
  hw  <- max(w/2*0.94, 0.30*w) * (if (big) 1.34 else 1.12)
  pad <- lh * (if (big) 0.60 else 0.34)
  yl  <- min(ycell - BH/2, ylab - lh*0.75) - pad
  yh  <- max(ycell + BH/2, ylab + lh*0.75) + pad
  th  <- seq(0, 2*pi, length.out = 181)
  polygon(xc + hw*cos(th), (yl+yh)/2 + ((yh-yl)/2)*sin(th),
          border = "black", lwd = lwd, col = NA)
}
kof <- function(key) nchar(sub("^[^ ]+ ", "", key))   # mechanisms in a cell key
## Contenders (DIC only; RECOVERY_jm.md, 3 Nov): a model not separable from its model type's best carries "=".
## Directional thresholds from the recovery studies, rounded (Andrew, 3 Nov): a model with FEWER mechanisms than the
## best is a contender within ADD_THR (a more complex wrong winner beat the truth by at most 40, studies 1-4, and 28,
## study 5); any other model (more mechanisms or the same number) within DROP_THR (the floor, 149, rounded). The rule
## is stated in README_display.md, not on the figure.
ADD_THR <- 50; DROP_THR <- 150
## Recovery study 5 rejects these one-mechanism models although they are within ADD_THR of the best: in their cells,
## one-mechanism data always put the one-mechanism model clearly ahead (rec5/recovery_rec5_cells.csv). Keys "task arch mech".
REC5_OUT <- c("Simon RDM b", "Simon LBA b")
## "=" as two white bars with a thin dark outline, legible on every fill (Andrew, 3 Nov); on a hatched (partially converged)
## cell a patch of the fill clears the lines behind it. cex is unused (kept for the call sites).
## sx, sy scale the bars (legend: larger, so the thin outline does not swallow a bar at screen zoom; Andrew, 4 Nov)
## Cells too (Andrew, 4 Nov: at 1x the outline filled whole bars black at screen zoom): sx 1.5, sy 1.8 by default.
eqmark <- function(xc, y, cex, fill = NULL, hatched = FALSE, sx = 1.5, sy = 1.8) {
  if (hatched) rect(xc - 0.0055*sx, y - min(BH*0.36*sy, BH*0.47), xc + 0.0055*sx, y + min(BH*0.36*sy, BH*0.47), col = fill, border = NA)
  for (o in c(-1, 1)) rect(xc - 0.0040*sx, y + o*BH*0.15*sy - BH*0.065*sy, xc + 0.0040*sx, y + o*BH*0.15*sy + BH*0.065*sy,
                           col = "white", border = "grey10", lwd = 0.5) }
NCONT <- NULL
## collect a tree's cells, then box its own minimum once the tree is drawn
TB <- NULL
tb_add <- function(arch, set, xc, ycell, ylab, w)
  TB[[length(TB)+1]] <<- list(v = getv(arch, set), k = paste(arch, paste(sort(set), collapse="")),
                              x = xc, yc = ycell, yl = ylab, w = w)
tb_draw <- function(LH) {
  vv <- sapply(TB, function(z) if (is.null(z$v) || is.na(z$v) || nobest(z$k)) Inf else z$v)
  if (all(!is.finite(vv))) { TB <<- list(); return(invisible()) }
  z <- TB[[which.min(vv)]]
  if (!identical(z$k, BEST)) bestbox(z$x, z$yc, z$yl, z$w, LH, lwd = 1)
  ## delta / CAF (Andrew, 4 Nov): "=" within the difference floor PNF of the type's best misfit, no direction
  { nc <- 0
    for (e in TB) { if (is.null(e$v) || is.na(e$v) || identical(e$k, z$k) || (MEASURE == "dic" && paste(TK, e$k) %in% REC5_OUT)) next
      thr <- if (MEASURE == "dic") (if (kof(e$k) < kof(z$k)) ADD_THR else DROP_THR) else PNF
      if (!is.na(thr) && e$v - z$v <= thr) { eqmark(e$x, e$yc, CXP, col_of(e$v), paste(TK, e$k) %in% RIDGE); nc <- nc + 1 } }
    NCONT[[sub(" .*$", "", z$k)]] <<- nc }
  i2 <- which(vapply(TB, function(z) kof(z$k) == 2, TRUE) & is.finite(vv))
  if (length(i2)) { z2 <- TB[[ i2[which.min(vv[i2])] ]]
    bestcircle(z2$x, z2$yc, z2$yl, z2$w, LH,
               lwd = if (!is.na(BEST2) && identical(z2$k, BEST2)) 3.8 else 1) }
  TB <<- list()
}
balance <- function(kids, parents) {
  cnt <- integer(length(parents)); own <- integer(length(kids))
  for (c in seq_along(kids)) { ok <- which(sapply(parents, function(p) all(kids[[c]] %in% p)))
    own[c] <- ok[which.min(cnt[ok])]; cnt[own[c]] <- cnt[own[c]] + 1 }
  own
}
lab4 <- function(z, a) sub(paste0("^", a, "_"), "", encode(paste(z, collapse=""), a))
labR <- function(z)    sub("^RDM_", "", encode(paste(z, collapse=""), "RDM"))

## ---- four-mechanism panel: 1 - 4 - 6 - 4 -----------------------------------
build4 <- function(a) {
  L <- lapply(4:1, function(k) combn(FOUR[[a]], k, simplify = FALSE)); own <- vector("list", 4)
  for (i in 1:3) { o <- balance(L[[i+1]], L[[i]]); ord <- order(o, seq_along(o))
                   L[[i+1]] <- L[[i+1]][ord]; own[[i+1]] <- o[ord] }
  list(L = L, own = own, lb = lapply(L, function(lv) sapply(lv, lab4, a = a)))
}
draw4 <- function(B, a, ytop, ybot, xl, xr, CX) {
  L <- B$L; own <- B$own; lb <- B$lb
  LH <- strheight("Xg", cex = CX); LOFF <- BH/2 + LH/2 + 0.008; LTOP <- LOFF + LH/2 + 0.005
  YL <- ytop - (0:3)*(ytop - ybot)/3; W <- (xr - xl)/6
  ctr <- function(n, sp) (xl+xr)/2 - sp/2 + (seq_len(n) - 0.5)*(sp/n)
  X <- list((xl+xr)/2, ctr(4, 6*W), ctr(6, 6*W), NULL)
  tgt <- seq(min(X[[3]]), max(X[[3]]), length.out = 4); used <- integer(0)
  X[[4]] <- numeric(4); own[[4]] <- integer(4)
  for (j in 1:4) { ok <- which(sapply(L[[3]], function(p) all(L[[4]][[j]] %in% p)) &
                               !(seq_along(L[[3]]) %in% used))
    p <- ok[which.min(abs(X[[3]][ok] - tgt[j]))]; used <- c(used, p)
    own[[4]][j] <- p; X[[4]][j] <- X[[3]][p] }
  for (i in 1:3) for (j in seq_along(L[[i+1]]))
    segments(X[[i]][own[[i+1]][j]], YL[i]-BH/2, X[[i+1]][j], YL[i+1]+LTOP, lwd=0.6, col="grey40")
  for (i in 1:4) for (j in seq_along(L[[i]])) {
    cell(X[[i]][j], YL[i], W, getv(a, L[[i]][[j]]), getst(a, L[[i]][[j]]))
    text(X[[i]][j], YL[i]+LOFF, lb[[i]][j], cex=CX)
    tb_add(a, L[[i]][[j]], X[[i]][j], YL[i], YL[i]+LOFF, W)
    if (isbest(a, L[[i]][[j]])) bestbox(X[[i]][j], YL[i], YL[i]+LOFF, W, LH) }
  tb_draw(LH)
  text((xl+xr)/2, YL[1] + LOFF + LH*1.9, archlab(a), cex = CX*1.5, font = 2)
}
## ---- race panel: three trees, 1 - 4 - 5 - (2 shared) -----------------------
drawR <- function(ytop, ybot, xl, xr, header, CX) {
  ARCHNOW <- header
  TB <<- list()
  LH <- strheight("Xg", cex = CX); LOFF <- BH/2 + LH/2 + 0.008; LTOP <- LOFF + LH/2 + 0.005
  YL <- ytop - (0:2)*(ytop - ybot)/3; YBOT <- ybot
  tw <- (xr - xl)/3; bottom <- list()
  for (ti in seq_along(ROOTS)) {
    root <- ROOTS[[ti]]; x0 <- xl + (ti-1)*tw
    L2 <- combn(root,3,simplify=FALSE); L3 <- combn(root,2,simplify=FALSE)
    drop <- which(sapply(L3, labR) %in% setdiff(c("Vci","Sci","Bci"), KEEP[ti]))
    kids <- L3[-drop]; own <- balance(kids, L2)
    o <- order(own, seq_along(own)); slots <- seq_along(L3)[-drop][o]
    for (sw in SWAP) { i1 <- match(sw[1], sapply(L3[slots], labR)); i2 <- match(sw[2], sapply(L3[slots], labR))
      if (!is.na(i1) && !is.na(i2)) { t<-slots[i1];slots[i1]<-slots[i2];slots[i2]<-t
                                      t<-o[i1];o[i1]<-o[i2];o[i2]<-t } }
    W <- (tw - 0.006)/5; xc <- x0 + tw/2; span <- 5*W
    x3 <- xc - span/2 + (seq_along(slots) - 0.5)*W
    x2 <- xc - span/2 + (seq_along(L2) - 0.5)*(span/length(L2))
    dl <- sapply(L3[slots], labR)
    haslink <- dl %in% BPAR
    below <- if (!any(haslink[seq(2,length(dl),2)])) seq(2,length(dl),2) else seq(1,length(dl),2)
    stopifnot(!any(haslink[below]))
    segments(xc, YL[1]-BH/2, x2, YL[2]+LTOP, lwd=0.6, col="grey40")
    ## a link into a cell whose label sits BELOW has no label to clear, so it can
    ## run all the way down to the bar
    for (j in seq_along(slots))
      segments(x2[own[o][j]], YL[2]-BH/2, x3[j],
               YL[3] + if (j %in% below) BH/2 + 0.003 else LTOP, lwd=0.6, col="grey40")
    cell(xc, YL[1], W, getv(ARCHNOW, root), getst(ARCHNOW, root))
    text(xc, YL[1]+LOFF, labR(root), cex=CX)
    tb_add(ARCHNOW, root, xc, YL[1], YL[1]+LOFF, W)
    if (isbest(ARCHNOW, root)) bestbox(xc, YL[1], YL[1]+LOFF, W, LH)
    for (j in seq_along(L2)) { cell(x2[j], YL[2], W, getv(ARCHNOW, L2[[j]]), getst(ARCHNOW, L2[[j]]))
      text(x2[j], YL[2]+LOFF, labR(L2[[j]]), cex=CX)
      tb_add(ARCHNOW, L2[[j]], x2[j], YL[2], YL[2]+LOFF, W)
      if (isbest(ARCHNOW, L2[[j]])) bestbox(x2[j], YL[2], YL[2]+LOFF, W, LH) }
    for (j in seq_along(slots)) cell(x3[j], YL[3], W, getv(ARCHNOW, L3[[slots[j]]]), getst(ARCHNOW, L3[[slots[j]]]))
    for (j in seq_along(slots)) {
      yl <- YL[3] + if (j %in% below) -LOFF else LOFF
      text(x3[j], yl, dl[j], cex=CX)
      tb_add(ARCHNOW, L3[[slots[j]]], x3[j], YL[3], yl, W)
      if (isbest(ARCHNOW, L3[[slots[j]]])) bestbox(x3[j], YL[3], yl, W, LH) }
    for (m in names(BPAR)) { j <- match(BPAR[m], dl)
      if (!is.na(j)) bottom[[m]] <- c(x3[j], W) }
  }
  for (m in c("V","v","S","s","B","b")) {
    x <- bottom[[m]][1]
    segments(x, YL[3]-BH/2, x, YBOT+LTOP, lwd=0.6, col="grey40")
    cell(x, YBOT, bottom[[m]][2], getv(ARCHNOW, m), getst(ARCHNOW, m))
    text(x, YBOT+LOFF, labR(m), cex=CX)
    tb_add(ARCHNOW, m, x, YBOT, YBOT+LOFF, bottom[[m]][2])
    if (isbest(ARCHNOW, m)) bestbox(x, YBOT, YBOT+LOFF, bottom[[m]][2], LH)
  }
  tb_draw(LH)                                   # one thin box for the whole panel
  text((xl+xr)/2, YL[1] + LOFF + strheight("Xg", cex=CX)*1.9, archlab(header),
       cex = CX*1.5, font = 2)
}
for (TK in TASKS) {
## one set of PDFs with the simplest names (Andrew, 2026-10-21): grid_<measure>_<task>.pdf; the paper's delta grid
## is the median one (PP_AGG=median; CAF is identical under either aggregate: accuracy is always a mean)
OUT <- file.path(OUTDIR, sprintf("grid_%s_%s.pdf", MEASURE, TK))
pdf(OUT, width = 19, height = 9.5); par(mar = c(0.3,0.3,0.3,0.3))
VAL <- VL[[TK]]; STATES <- character(0)
## Limits come from the cells actually DRAWN on the page, not from every
## value in the ensemble: the overall best is sometimes a large-k model
## that has no cell (Simon DIC's minimum of 0 is one), which left the best
## drawn cell at ramp 37 instead of 0. Per page, drawn-only: the best cell
## is always darkest blue and the worst always darkest red.
## For DIC the values are differences from a reference, so re-reference them to
## the best cell actually on the page: otherwise the zero belongs to a large-k
## model the tree never shows and nothing reads 0.
if (MEASURE %in% c("dic","jdic")) {
  dmin <- min(VAL[names(VAL) %in% DRAWN & !nobest(names(VAL))], na.rm = TRUE); VAL <- VAL - dmin }
{ vd <- VAL[names(VAL) %in% DRAWN]; vb <- vd[!nobest(names(vd))]; BEST <- names(vb)[which.min(vb)]
  v2 <- vb[vapply(names(vb), kof, 1L) == 2]
  BEST2 <- if (length(v2)) names(v2)[which.min(v2)] else NA_character_
  LO <- min(vb); HI <- max(vd); PNF <- nfloor(MEASURE, TK)
  if (!(MEASURE %in% c("delta","caf"))) {
    nf <- nfloor(MEASURE, TK); LOGOFF <- if (is.na(nf)) 1 else max(1, nf) } }
plot(NA, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "")
XL <- 0.012; UW <- 1 - 2*XL; CW <- UW/4          # four equal columns
B4 <- lapply(names(FOUR), build4); names(B4) <- names(FOUR)
## one label size: tightest of the top row (pitch W and 1.5W) and bottom (W and 2W)
Wt <- (CW - 0.010)/6; Wb <- (2*CW/3 - 0.006)/5
CX <- MAXCEX
for (a in names(FOUR)) CX <- min(CX, 0.90*Wt/max(strwidth(B4[[a]]$lb[[3]], cex=1)),
                                     0.90*1.5*Wt/max(strwidth(B4[[a]]$lb[[2]], cex=1)))
CX <- min(CX, 0.90*2*Wb/strwidth("BciVc", cex=1), 0.90*1.25*Wb/strwidth("VciSc", cex=1))
cat(sprintf("top cell %.4f, bottom cell %.4f, cex %.2f (cap %.2f)\n", Wt, Wb, CX, MAXCEX)); CXP <- CX
for (i in seq_along(FOUR)) { a <- names(FOUR)[i]
  draw4(B4[[a]], a, 0.93, 0.60, XL + (i-1)*CW + 0.005, XL + i*CW - 0.005, CX) }
drawR(0.43, 0.06, XL,          XL + 2*CW - 0.01, "RDM", CX)
drawR(0.43, 0.06, XL + 2*CW + 0.01, XL + 4*CW,   "LBA", CX)
## legend in the free top-right cell
lx <- XL + 3*CW + CW*0.10; ly <- seq(0.62, 0.875, length.out = 101)
for (k in 1:100) rect(lx, ly[k], lx + 0.016, ly[k+1], col = RAMP[k], border = NA)
rect(lx, ly[1], lx + 0.016, ly[101], border = "grey40", lwd = 0.7)
tk <- switch(MEASURE, delta = c(7,10,20,40,70,110), caf = c(0.7,1,2,4,8,17),
             c(0, 300, 1000, 3000, 10000))
for (v in tk) {
  yy <- ly[1] + (tf(v) - tf(LO))/(tf(HI) - tf(LO)) * (ly[101]-ly[1])
  if (yy >= ly[1] && yy <= ly[101]) {
    segments(lx + 0.016, yy, lx + 0.019, yy, lwd = 0.7, col = "grey40")
    text(lx + 0.022, yy, v, adj = 0, cex = CX*0.85) }
}

text(lx - 0.004, 0.945, switch(MEASURE, delta = "delta misfit", caf = "CAF misfit",
     dic = "conditional DIC", jdic = "joint DIC"), cex = CX*1.0, font = 2, adj = 0)
text(lx - 0.004, 0.915, switch(MEASURE, delta = "% of effect", caf = "percentage pts",
     "diff from best"), cex = CX*0.85, adj = 0)
## Variant legend, to the right of the thermometer. Labels are split over two
## lines so the block stays narrow.
vx <- lx + 0.085; bw <- 0.019; bh <- 0.012; LC <- CX*0.80
text(vx - 0.004, 0.945, "Variant", cex = CX*1.0, font = 2, adj = 0)
## Entries are listed only when the page has such a cell (3 Nov: no fit is running now, and "not applicable" cells
## exist only on the Stroop page), stacked from the top.
yy <- 0.880; nxt <- function() { y <- yy; yy <<- yy - 0.050; y }
lab2 <- function(y, a, b) { text(vx+bw+0.005, y+0.010, a, adj = 0, cex = LC); text(vx+bw+0.005, y-0.010, b, adj = 0, cex = LC) }
## order (Andrew, 3 Nov): best, not separable, best 2 mechanism, then partially converged, did not converge, [not applicable, fit running]
y <- nxt(); rect(vx, y-bh, vx+bw, y+bh, col = NA, border = "black", lwd = 3.0); text(vx+bw+0.005, y, "best", adj = 0, cex = LC)
{ y <- nxt(); rect(vx, y-bh*1.15, vx+bw, y+bh*1.15, col = "grey62", border = NA); eqmark(vx+bw/2, y, CX, sx = 1.6, sy = 1.9)
  lab2(y, "not separable", "from best") }
y <- nxt(); local({ th <- seq(0, 2*pi, length.out = 181)
  polygon(vx+bw/2 + (bw/2)*1.15*cos(th), y + bh*1.15*sin(th), border = "black", lwd = 3.0, col = NA) })
lab2(y, "best 2", "mechanism")
yy <- yy - 0.010                                          # a small gap between the markers and the cell states
if ("rg" %in% STATES) { y <- nxt(); rect(vx, y-bh, vx+bw, y+bh, col = "grey92", border = NA); ridgehatch(vx, y-bh, vx+bw, y+bh)
  text(vx+bw+0.005, y, "partially converged", adj = 0, cex = LC) }
if ("nc" %in% STATES) { y <- nxt(); rect(vx, y-bh, vx+bw, y+bh, col = "grey82", border = "grey62", lwd = 0.4); lab2(y, "did not", "converge") }
if ("na" %in% STATES) { y <- nxt(); rect(vx, y-bh, vx+bw, y+bh, col = NA, border = "grey55", lwd = 0.5)
  segments(vx, y-bh, vx+bw, y+bh, col = "grey55", lwd = 0.5); segments(vx, y+bh, vx+bw, y-bh, col = "grey55", lwd = 0.5)
  lab2(y, "not", "applicable") }
if ("rn" %in% STATES) { y <- nxt(); rect(vx, y-bh, vx+bw, y+bh, col = "grey95", border = "grey50", lwd = 0.5, lty = 2); lab2(y, "fit", "running") }
## Smallest difference worth reading (delta and CAF): the misfit a perfect model would still show, from resampling
## the subjects. DIC uses the "=" contender rule instead (README_display.md).
NF <- nfloor(MEASURE, TK)
if (MEASURE != "dic" && !is.na(NF)) {
  text(vx, y - 0.060, "Difference floor:", adj = 0, cex = LC, font = 2)
  text(vx, y - 0.085, sprintf("~%s", fmtv(NF)), adj = 0, cex = LC, font = 2)
}
{ cat(TK, MEASURE, "contenders per model type:", paste(names(NCONT), unlist(NCONT), collapse = ", "), "\n"); NCONT <- NULL }
invisible(dev.off()); cat("wrote", OUT, "\n")
}
if (any(NFCELL)) cat("cells with no joint-move fit (drawn 'did not converge'):", names(NFCELL)[NFCELL], "\n")
