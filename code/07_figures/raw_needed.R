## The fits whose simulated trials the EMC2-drawn delta / CAF figures plot, and their pp_raw files:
##   Rscript code/07_figures/raw_needed.R [cores]   -> derived/tables/figure_fits.csv; runs 06_predict/pp_subj.R raw for
##                                                     any label without derived/pp_raw/<label>_subj.rds
## best5: plot_best5.R's own picks (LIST_ONLY=1); models6: RDM VciSci, ViSi, BciSci, BiSi and LNR MciSci, MiSi per task;
## models6 mix row 3 (the paper's Figure 11): LNR Mi, RDM Vi and RDM Bi for Stroop
## and ANT; choice: the five full simple-choice variants (their pp_raw is made by plot_choice.R sim). The Figure 12
## knockouts are simulated from the RDM ViSi fits (06_predict/pp_null.R, derived/pp_null/).
## figure_fits.csv is also the list of fits published with the online data record.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
cores <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "4"
D <- read_tab("display_jm.csv")
lab_of <- function(spec) D$label[match(spec, D$spec)]
b5 <- system2("Rscript", shQuote(file.path(CODE, "07_figures", "plot_best5.R")), stdout = TRUE, env = "LIST_ONLY=1")
b5 <- unique(sub("^NEED ([^ ]+).*", "\\1", grep("^NEED ", b5, value = TRUE)))
m6 <- lab_of(as.vector(outer(TASKS, c("RDM_VciSci", "RDM_ViSi", "RDM_BciSci", "RDM_BiSi", "LNR_MciSci", "LNR_MiSi"), paste, sep = "_")))
mx <- lab_of(as.vector(outer(c("Stroop", "ANT"), c("LNR_Mi", "RDM_Vi", "RDM_Bi"), paste, sep = "_")))
ch <- unname(CHOICE_FULL)
F <- rbind(data.frame(label = b5, figure = "best5"), data.frame(label = m6, figure = "models6"),
           data.frame(label = mx, figure = "models6mix"),
           data.frame(label = ch, figure = "choice"))
F <- aggregate(figure ~ label, F, function(x) paste(sort(unique(x)), collapse = " "))
F$pp_raw <- file.exists(pp_raw_path(F$label))
write_out(F, tab_path("figure_fits.csv")); print(F, row.names = FALSE)
need <- F$label[!F$pp_raw & !grepl("^Choice_", F$label)]
if (length(need)) system2("Rscript", c(shQuote(file.path(CODE, "06_predict", "pp_subj.R")), "raw", cores, need))
if (any(!F$pp_raw & grepl("^Choice_", F$label))) system2("Rscript", c(shQuote(file.path(CODE, "07_figures", "plot_choice.R")), "sim", cores))
