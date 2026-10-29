options(stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
stage <- if(length(args)) args[1] else "check"
f <- sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
code_dir <- dirname(normalizePath(f))
root <- dirname(code_dir)
data_dir <- Sys.getenv("CEBP_DATA_DIR",file.path(root,"data"))
results_dir <- Sys.getenv("CEBP_RESULTS_DIR",file.path(root,"results"))
dir.create(results_dir,recursive=TRUE,showWarnings=FALSE)
stages <- c(recode="clinical/recode.R",baseline="clinical/baseline.R",screening="clinical/screening.R",multivariable="clinical/multivariable.R",subgroups="clinical/subgroups.R",screen_plot="clinical/screen_plot.R",ukb_score="ukb/score.R",ukb_dimensions="ukb/dimensions.R")
if(stage=="check") {
 files <- list.files(code_dir,pattern="[.]R$",recursive=TRUE,full.names=TRUE)
 invisible(lapply(files,parse))
 checks <- data.frame(input=c("clinical_survival.csv","ukb_analysis.csv"))
 checks$available <- file.exists(file.path(data_dir,checks$input))
 write.csv(checks,file.path(results_dir,"input_check.csv"),row.names=FALSE)
 capture.output(sessionInfo(),file=file.path(results_dir,"session_info.txt"))
 cat("Syntax check passed. This is not a full reproducibility run.
")
 print(checks)
 quit(status=0)
}
if(stage %in% c("clinical","recode","baseline","screening","multivariable","subgroups","screen_plot")) {
 source(file.path(code_dir,"clinical","engine.R"))
 run_clinical(stage,code_dir,data_dir,results_dir)
 quit(status=0)
}
if(!stage %in% names(stages))stop("Unknown stage: ",stage)
required <- "ukb_analysis.csv"
if(!file.exists(file.path(data_dir,required)))stop("Provide authorized input: ",required)
stage_dir <- file.path(results_dir,stage)
dir.create(stage_dir,recursive=TRUE,showWarnings=FALSE)
warning_log <- character()
withCallingHandlers(source(file.path(code_dir,stages[[stage]]),local=environment()), warning=function(w){warning_log <<- c(warning_log,conditionMessage(w))})
write.csv(data.frame(warning=warning_log),file.path(stage_dir,"warnings.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(stage_dir,"session_info.txt"))
