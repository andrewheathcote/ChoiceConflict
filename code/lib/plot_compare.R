## plot_compare(): the information-criterion panel of the simple-choice figure, from NCI/Functions.R (lines 1-37).
plot_compare <- function(ict,do=c("WDM","DDM","RDM","LNR","LBA"),ylim=NULL,
  use="DIC",lpos="top",spaceby=.1,mcex=.75,
  tnames=c(f="Bias",u="Unbiased",s="Bias+Eq.SD",us="Unbiased+Eq.SD")) {

  parseNames <- function(nams) {
    out <- substr(nams,1,3)
    attr(out,"type") <- substring(nams,4)
    out
  }

  nams <- parseNames(rownames(ict))
  type <- attr(nams,"type")
  isin <- nams %in% do
  ict <- ict[isin,]
  ic <- ict[,use]
  icnams <- nams[isin]; type <- type[isin]
  odr <- order(setNames(c(1:length(do)),do)[nams])
  nams <- nams[odr]; type <- type[odr]; ic <- ic[odr]
  n <- table(nams)[do]
  space <- rep(0,sum(n))
  space[1+cumsum(n[-length(n)])] <- spaceby
  types <- unique(type)
  grays <- (as.numeric(factor(type,levels=types)))^2
  grays <- gray(1-grays/max(grays))
  ic <- ic-min(ic)
  if (is.null(ylim)) ylim <- c(0,max(ic))
  ylab <- paste0(use,"- min(",use,")")
  mid <- barplot(ic,space=space,col=grays,ylim=ylim,ylab=ylab)
  mtext(do,1,line=1,cex=mcex,at=tapply(mid,factor(nams,levels=do),mean))
  if(!is.null(tnames)) {
    types[1] <- "f"
    tnames <- tnames[types]
    legend(lpos,tnames,lty=rep(1,length(types)),lwd=6,
         col=unique(grays),bty="n")
  }
}

