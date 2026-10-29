# Integration test with public survival::lung data; not study participant data.
f <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(dirname(normalizePath(f)))
work <- Sys.getenv("CEBP_TEST_DIR",tempfile("clinical_test_"));dir.create(work,recursive=TRUE,showWarnings=FALSE)
x <- survival::lung
x <- data.frame(id=sprintf("p%04d",seq_len(nrow(x))),time=x$time,status=as.integer(x$status==2),age=x$age,sex=x$sex,ph.ecog=x$ph.ecog)
write.csv(x,file.path(work,"clinical_survival.csv"),row.names=FALSE)
clinical_cox_config <- list(input="clinical_survival.csv",time_unit="days",time_origin="Public lung dataset origin",endpoint_definition="Death in public lung test data",factors=list(sex=c(1,2)),screening=c("age","sex"),adjustment=c("age","sex"),exposures=c("age","sex"),subgroups="sex",baseline=c("age","sex"))
configfile<-file.path(work,"config.R")
writeLines(c("clinical_cox_config <-",capture.output(dput(clinical_cox_config))),configfile)
Sys.setenv(CEBP_DATA_DIR=work,CEBP_RESULTS_DIR=file.path(work,"results"),CEBP_CLINICAL_CONFIG=configfile)
runner<-file.path(root,"code","main.R")
rscript<-file.path(R.home("bin"),"Rscript.exe")
if(!file.exists(rscript))rscript<-file.path(R.home("bin"),"Rscript")
run <- function(stage)system2(rscript,c(shQuote(runner),stage),stdout=TRUE,stderr=TRUE)
log<-run("clinical");cat(paste(log,collapse="\n"),"\n")
stopifnot(is.null(attr(log,"status")))
out<-file.path(work,"results","clinical_cox")
z<-read.csv(file.path(out,"adjusted_HR.csv"))
x$sex<-factor(x$sex,levels=c(1,2))
fitted<-survival::coxph(survival::Surv(time,status==1)~age+sex,x,ties="efron")
expected<-exp(coef(fitted))
stopifnot(all(abs(z$HR-expected[z$term])<1e-10))
ci<-exp(confint(fitted));stopifnot(all(abs(z$lower-ci[z$term,1])<1e-8),all(abs(z$upper-ci[z$term,2])<1e-8))
stopifnot(file.info(file.path(out,"screening_HR.pdf"))$size>0)
ints<-read.csv(file.path(out,"interaction_tests.csv"));stopifnot(nrow(ints)==1,is.finite(ints$P))
for(s in c("check","recode","baseline","screening","multivariable","subgroups","screen_plot")) {
 log<-run(s);stopifnot(is.null(attr(log,"status")))
}
# Competing events are censored at their own time in cause-specific Cox.
x$status[seq(1,nrow(x),by=13)]<-2L
x$age[1]<-NA_real_
write.csv(x,file.path(work,"clinical_survival.csv"),row.names=FALSE)
log<-run("multivariable");stopifnot(is.null(attr(log,"status")))
z<-read.csv(file.path(out,"adjusted_HR.csv"))
fitted<-survival::coxph(survival::Surv(time,status==1)~age+sex,x,ties="efron")
stopifnot(all(abs(z$HR-exp(coef(fitted))[z$term])<1e-10),all(z$omitted_missing==1))
x$time[1]<-NA_real_;write.csv(x,file.path(work,"clinical_survival.csv"),row.names=FALSE)
log<-suppressWarnings(run("multivariable"));stopifnot(!is.null(attr(log,"status")))
x$time[1]<-survival::lung$time[1]
x$id[2]<-x$id[1];write.csv(x,file.path(work,"clinical_survival.csv"),row.names=FALSE)
log<-suppressWarnings(run("multivariable"));stopifnot(!is.null(attr(log,"status")))
cat("PASS: complete runner, stage aliases, direct Cox HR/CI agreement, subgroup interaction, forest plot, competing-event coding, missing covariate counts, missing-time and duplicate-ID rejection.\n")
cat("Test output:",work,"\n")

