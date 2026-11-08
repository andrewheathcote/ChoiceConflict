## Case-safe canonical names, grouped by parameter.
## B = threshold across ALL architectures (the DDM/WDM boundary `a` included).
## The diffusion start point `z` is a threshold asymmetry, so it is written as
## the i contrast on B rather than as a separate parameter: a -> Bc, z -> Bi,
## and a model with both -> Bci. That makes WDM/DDM directly comparable with the
## race models.
##   <PARAM><contrasts> where contrast c = congruency, i = iMd.
## A parameter carrying both contrasts is written once: VvSs -> VciSci.
## Parameter order is fixed (threshold, rate, noise), so names are canonical.
PARAM <- list(                       # mechanism letter -> (parameter, contrast)
  RDM = list(B=c("B","c"), b=c("B","i"), V=c("V","c"), v=c("V","i"), S=c("S","c"), s=c("S","i")),
  LBA = list(B=c("B","c"), b=c("B","i"), V=c("V","c"), v=c("V","i"), S=c("S","c"), s=c("S","i")),
  LNR = list(M=c("M","c"), m=c("M","i"), S=c("S","c"), s=c("S","i")),
  WDM = list(a=c("B","c"), z=c("B","i"), v=c("V","c"), s=c("S","c"), V=c("SV","c")),
  DDM = list(a=c("B","c"), z=c("B","i"), v=c("V","c"), s=c("S","c"), V=c("SV","c")))
PORD <- list(RDM=c("B","V","S"), LBA=c("B","V","S"), LNR=c("M","S"),
             WDM=c("B","V","S","SV"), DDM=c("B","V","S","SV"))
encode <- function(mech, arch) {
  if (is.na(mech) || !nzchar(mech)) return(NA_character_)
  ls <- unique(strsplit(mech, "")[[1]])
  tab <- PARAM[[arch]]
  if (!all(ls %in% names(tab))) return(NA_character_)
  pc <- do.call(rbind, lapply(ls, function(l) tab[[l]]))
  out <- ""
  for (p in PORD[[arch]]) {
    cs <- pc[pc[,1] == p, 2]
    if (length(cs)) out <- paste0(out, p, paste(c("c","i")[c("c","i") %in% cs], collapse=""))
  }
  paste0(arch, "_", out)
}
