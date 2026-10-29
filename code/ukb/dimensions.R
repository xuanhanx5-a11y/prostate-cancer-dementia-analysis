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
  grip = bad(Grip_h, "Abnormal"),
  fev1 = bad(fev1_best_l_cat, "Abnormal"),
  fvc = bad(fvc_best_l_cat, "Abnormal"),
  obesity = bad(BMI3_h, "28 or higher"),
  diabetes = as.integer(Glucose3_h == "Abnormal"),
  hypertension = bad(Hypertension_h, "Yes"),
  glucose = bad(glucose_cat, c("Borderline", "Abnormal")),
  triglycerides = bad(triglycerides_cat, c("Borderline", "Abnormal")),
  tyg = bad(tyg_index_component_cat, "High"),
  creatinine = bad(creatinine_cat, "Abnormal"),
  haemoglobin = bad(hemoglobin_cat, "Abnormal"),
  crp = bad(crp_cat, "High"),
  wbc = bad(wbc_cat, "Abnormal"),
  monocyte = bad(monocyte_cat, "Abnormal"),
  ast = bad(ast_cat, c("Borderline", "Abnormal")),
  ggt = bad(ggt_cat, c("Borderline", "Abnormal"))
)]
ukb[, patient_functional := rowMeans(.SD),
  .SDcols = c("smoking", "sedentary", "sleep", "falls", "grip", "fev1", "fvc")]
ukb[, clinical_biological := rowMeans(.SD), .SDcols = c(
  "obesity", "diabetes", "hypertension", "glucose", "triglycerides", "tyg",
  "creatinine", "haemoglobin", "crp", "wbc", "monocyte", "ast", "ggt"
)]
ukb[, `:=`(
  domain_behavioural = rowMeans(.SD) >= 2 / 3,
  domain_physical = FALSE,
  domain_cardiometabolic = FALSE,
  domain_systemic = FALSE
), .SDcols = c("smoking", "sedentary", "sleep")]
ukb[, domain_physical := rowMeans(.SD) >= 1 / 4,
  .SDcols = c("falls", "grip", "fev1", "fvc")]
ukb[, domain_cardiometabolic := rowMeans(.SD) >= 2 / 6,
  .SDcols = c("obesity", "diabetes", "hypertension", "glucose", "triglycerides", "tyg")]
ukb[, domain_systemic := rowMeans(.SD) >= 2 / 7,
  .SDcols = c("creatinine", "haemoglobin", "crp", "wbc", "monocyte", "ast", "ggt")]


ukb[, `:=`(pf_z=zscore(patient_functional),cb_z=zscore(clinical_biological),education_f=factor(Education2_h),ses_f=factor(SES_baseline_cat),rp_f=factor(rp_any))]
fit <- coxph(Surv(time,event) ~ pf_z * cb_z + age_i0 + education_f + ses_f + rp_f,data=ukb,ties="efron",x=TRUE)
z <- summary(fit)
fwrite(data.table(term=rownames(z$coefficients),HR=z$conf.int[,1],lower=z$conf.int[,3],upper=z$conf.int[,4],P=z$coefficients[,5]),file.path(stage_dir,"joint_dimensions.csv"))
r1 <- coxph(Surv(time,event) ~ cb_z + age_i0 + education_f + ses_f + rp_f,data=ukb,ties="efron")
r2 <- coxph(Surv(time,event) ~ pf_z + age_i0 + education_f + ses_f + rp_f,data=ukb,ties="efron")
out <- rbindlist(lapply(seq_along(list(r1,r2)),function(i){r<-list(r1,r2)[[i]];lr<-2*as.numeric(logLik(fit)-logLik(r));df<-attr(logLik(fit),"df")-attr(logLik(r),"df");data.table(removed=c("Patient/functional","Clinical/biological")[i],LR=lr,DF=df,P=pchisq(lr,df,lower.tail=FALSE),Delta_AIC=AIC(r)-AIC(fit))}))
fwrite(out,file.path(stage_dir,"dimension_ablation.csv"))
