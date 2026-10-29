suppressPackageStartupMessages({library(data.table); library(survival)})
ukb <- fread(file.path(data_dir, "ukb_analysis.csv"))
bad <- function(x, values) as.integer(as.character(x) %chin% values)
zscore <- function(x) as.numeric(scale(x))
cl_beh <- c("smoking", "sedentary", "sleep", "falls")
cl_bio <- c("obesity", "diabetes", "hypertension", "glucose", "creatinine", "haemoglobin", "crp", "triglycerides", "tyg", "wbc", "monocyte", "ast", "ggt")
ukb[, `:=`(
  event = as.integer(mh_dementia_ever == 1L),
  time = as.numeric(mh_dementia_delta_days_vs_index),
  smoking = bad(Smoking3_h, c("Former", "Current")),
  sedentary = bad(sedentary_cat, "High"),
  sleep = as.integer(!grepl("Normal", as.character(sleep_cat))),
  falls = as.integer(grepl("2 falls", as.character(Falls3_h), fixed = TRUE)),
  obesity = bad(BMI3_h, "28 or higher"),
  diabetes = as.integer(Glucose3_h == "Abnormal"),
  hypertension = bad(Hypertension_h, "Yes"),
  glucose = bad(glucose_cat, c("Borderline", "Abnormal")),
  creatinine = bad(creatinine_cat, "Abnormal"),
  haemoglobin = bad(hemoglobin_cat, "Abnormal"),
  crp = bad(crp_cat, "High"),
  triglycerides = bad(triglycerides_cat, c("Borderline", "Abnormal")),
  tyg = bad(tyg_index_component_cat, "High"),
  wbc = bad(wbc_cat, "Abnormal"),
  monocyte = bad(monocyte_cat, "Abnormal"),
  ast = bad(ast_cat, c("Borderline", "Abnormal")),
  ggt = bad(ggt_cat, c("Borderline", "Abnormal"))
)]
ukb[, core_score_locked := rowMeans(.SD), .SDcols = c(cl_beh, cl_bio)]
ukb[, core_z_locked := zscore(core_score_locked)]


fit <- coxph(Surv(time,event) ~ core_z_locked + age_i0 + factor(Education2_h) + factor(SES_baseline_cat) + factor(rp_any), data=ukb, ties="efron", x=TRUE)
z <- summary(fit)
fwrite(data.table(term=rownames(z$coefficients), HR=z$conf.int[,1], lower=z$conf.int[,3], upper=z$conf.int[,4], P=z$coefficients[,5]), file.path(stage_dir,"score_cox.csv"))
fwrite(data.table(N=fit$n,events=fit$nevent,C_index=z$concordance[1],C_SE=z$concordance[2],AIC=AIC(fit)),file.path(stage_dir,"model_summary.csv"))
