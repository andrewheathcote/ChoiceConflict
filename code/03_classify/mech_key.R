## One mechanism key for every displayed model: spec -> mechanism letters (name_scheme.R's codes, e.g. "Vs").
## Mechanisms are never taken from file names. Source: the jm index (frozen/index_jm.csv, by spec); the campaign DDMs, which
## are not in it, from the legacy archive's audit (frozen/audit_mech.csv, by task and variant).
##   Rscript code/03_classify/mech_key.R      -> derived/tables/mech_key.csv
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
D  <- read_tab("display_jm_base.csv")
IX <- read_tab("frozen/index_jm.csv"); AU <- read_tab("frozen/audit_mech.csv")
K <- data.frame(spec = D$spec, task = D$task, arch = D$arch, variant = sub("^[^_]+_", "", D$spec))
K$mech <- IX$mech[match(K$spec, IX$spec)]
k <- is.na(K$mech); K$mech[k] <- AU$mech[match(paste(K$task[k], K$variant[k]), paste(AU$task, AU$variant))]
K$source <- ifelse(k, "audit_ensembles (task, variant)", "index_jm (spec)")
if (anyNA(K$mech)) stop("no mechanism for: ", paste(K$spec[is.na(K$mech)], collapse = ", "))
write_out(K, tab_path("mech_key.csv")); print(table(K$source))
