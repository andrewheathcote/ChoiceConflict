## A readable table of every archived model's specification, from inputs/specs/<spec>.rds:
##   Rscript code/01_specs/describe_specs.R     -> derived/tables/specs_designs.csv
## One row per spec: formulas, contrasts, constants, order restrictions (pre-transforms), chains, number of subject
## parameters and their names, data file. A spec is everything a fit needs besides the sampler settings: 02_fit/fit.R
## builds the design from it with lib/transplant.R::build_emc(), which checks the parameter names and the prior.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
dep <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
cm <- function(m) paste(apply(m, 1, function(r) paste(format(r, digits = 3), collapse = ",")), collapse = "; ")
files <- list.files(file.path(INPUTS, "specs"), "[.]rds$", full.names = TRUE)
D <- do.call(rbind, lapply(files, function(f) { s <- readRDS(f)
  data.frame(spec = sub("[.]rds$", "", basename(f)), task = s$task, arch = s$arch, model = s$model,
             formulas = paste(vapply(s$flist, dep, ""), collapse = " | "),
             contrasts = paste(unlist(lapply(names(s$clist), function(p) { x <- s$clist[[p]]
               if (is.matrix(x)) sprintf("%s=[%s]", p, cm(x)) else paste(sprintf("%s:%s=[%s]", p, names(x), vapply(x, cm, "")), collapse = " ") })), collapse = " "),
             constants = if (length(s$consts)) paste(sprintf("%s=%s", names(s$consts), format(s$consts, digits = 4)), collapse = " ") else "",
             restrictions = if (length(s$pre$func)) paste(sprintf("%s=%s", names(s$pre$func), s$pre$func), collapse = " ") else "",
             chains = s$n_chains, n_par = length(s$stored_par_names), parameters = paste(s$stored_par_names, collapse = " "),
             data_file = basename(s$data_file)) }))
write_out(D, tab_path("specs_designs.csv")); cat(nrow(D), "specs described\n")
