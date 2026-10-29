suppressPackageStartupMessages(library(survival))
f<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]);code<-dirname(dirname(normalizePath(f)));root<-dirname(code)
source(file.path(code,"general/methods.R"));source(file.path(code,"general/config.R"))
data_dir<-Sys.getenv("CEBP_DATA_DIR",file.path(root,"data"));out<-file.path(Sys.getenv("CEBP_RESULTS_DIR",file.path(root,"results")),"general");dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-read.csv(file.path(data_dir,config$input),stringsAsFactors=FALSE);d<-check_survival(d)
for(v in config$factors)d[[v]]<-factor(d[[v]])
write_out<-function(x,name)write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE)
# Fail on incomplete analysis inputs; run multiple_impute explicitly rather than silently drop records.
if(length(config$indicators)){s<-score_from_indicators(d,config$indicators);d$score<-s$z;write_out(data.frame(center=s$center,scale=s$scale),"score_scaling")}
fit<-cox_model(d,config$cox_rhs);write_out(cox_table(fit),"dementia_cox");write_out(ph_test(fit),"proportional_hazards")
death<-cox_model(d,config$cox_rhs,2);write_out(cox_table(death),"death_cox")
fg<-fine_gray(d,config$fg_rhs);write_out(fg$table,"fine_gray")
if(length(config$screen))write_out(screen_cox(d,config$screen,config$screen_adjust),"screening")
if(length(config$period_cuts)){z<-period_cox(d,config$cox_rhs,config$period_cuts);write_out(do.call(rbind,lapply(seq_along(z),function(i)cbind(period=i,cox_table(z[[i]])))),"period_models")}
if(length(config$model_rhs)){z<-incremental_models(d,config$model_rhs);write_out(z$table,"model_comparison")}
for(lag in config$lags)write_out(cox_table(landmark_cox(d,config$cox_rhs,lag)),paste0("landmark_",lag))
if(!is.null(config$group)){aj<-aj_curve(d,config$group);saveRDS(aj,file.path(out,"aj_curves.rds"));pdf(file.path(out,"cumulative_incidence.pdf"));plot(aj,xlab=config$time_unit,ylab="State probability");dev.off()}
writeLines("Net risk, competing-risk absolute risk and bootstrap require an explicit reference population/horizon. Invoke the documented functions; no reference values are silently chosen.",file.path(out,"absolute_risk_setup.txt"))
capture.output(sessionInfo(),file=file.path(out,"session_info.txt"))
