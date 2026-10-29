run_clinical <- function(stage,code_dir,data_dir,results_dir) {
 source(file.path(code_dir,"general","methods.R"),local=TRUE)
 config_path <- Sys.getenv("CEBP_CLINICAL_CONFIG",file.path(code_dir,"clinical","cox_config.R"))
 source(config_path,local=TRUE)
 cfg <- clinical_cox_config
 for(v in c("time_unit","time_origin","endpoint_definition"))
   if(is.null(cfg[[v]]) || length(cfg[[v]])!=1L || !nzchar(cfg[[v]])) stop("Specify ",v," in clinical configuration")
 out <- file.path(results_dir,"clinical_cox")
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 d <- read.csv(file.path(data_dir,cfg$input),check.names=FALSE,colClasses=c(id="character"))
 check_survival(d)
 if(!is.numeric(d$time)||any(!is.finite(d$time)))stop("Require finite numeric time")
 if(any(!nzchar(d$id)))stop("Empty ID")
 needed <- unique(c(cfg$screening,cfg$adjustment,cfg$exposures,cfg$subgroups,cfg$baseline,names(cfg$factors)))
 if(!all(needed %in% names(d)))stop("Missing configured columns: ",paste(setdiff(needed,names(d)),collapse=", "))
 if(any(grepl("`",needed,fixed=TRUE)))stop("Backticks in variable names are unsupported")
 for(v in names(cfg$factors)) {
   lev <- cfg$factors[[v]]
   if(any(!is.na(d[[v]]) & !d[[v]] %in% lev))stop("Unexpected factor level: ",v)
   d[[v]] <- factor(d[[v]],levels=lev)
 }
 quote_name <- function(v) if(length(v)) paste0("`",v,"`") else character()
 warnings <- character()
 capture_fit <- function(expr) withCallingHandlers(expr,warning=function(w){warnings <<- c(warnings,conditionMessage(w))})
 fit_one <- function(v,adjust,input=d,label="adjusted") {
   columns <- unique(c("id","time","status",v,adjust))
   x <- input[complete.cases(input[,columns,drop=FALSE]),columns,drop=FALSE]
   if(!nrow(x)||sum(x$status==1)==0)stop("No usable dementia events: ",v)
   rhs <- paste(quote_name(c(v,setdiff(adjust,v))),collapse=" + ")
   before <- length(warnings)
   f <- capture_fit(cox_model(x,rhs))
   mm <- model.matrix(f)
   terms <- colnames(mm)[attr(mm,"assign")==1L]
   z <- cox_table(f); z <- z[z$term %in% terms,,drop=FALSE]
   z$exposure <- v; z$model <- label; z$omitted_missing <- nrow(input)-nrow(x)
   z$warnings <- paste(tail(warnings,length(warnings)-before),collapse="; ")
   z$finite_estimates <- is.finite(z$HR)&is.finite(z$lower)&is.finite(z$upper)
   list(table=z,fit=f)
 }
 save_table <- function(z,name)write.csv(z,file.path(out,paste0(name,".csv")),row.names=FALSE)
 selected <- function(s)stage %in% c("clinical",s)
 if(selected("recode"))save_table(d,"analysis_ready")
 if(selected("baseline")) {
   if(length(cfg$baseline))save_table(baseline_compare(transform(d,event_group=as.integer(status==1)),"event_group",cfg$baseline),"baseline")
   save_table(data.frame(N=nrow(d),dementia=sum(d$status==1),death_before_dementia=sum(d$status==2),censored=sum(d$status==0)),"endpoint_counts")
 }
 if(selected("screening")) {
   if(!length(cfg$screening))stop("Configure screening variables")
   fits <- lapply(cfg$screening,fit_one,adjust=character(),label="univariable")
   z <- do.call(rbind,lapply(fits,`[[`,"table")); z$FDR_contrast <- p.adjust(z$P,"BH")
   save_table(z,"screening_HR")
 }
 if(selected("multivariable")) {
   if(!length(cfg$exposures))stop("Configure exposures")
   fits <- lapply(cfg$exposures,fit_one,adjust=cfg$adjustment)
   save_table(do.call(rbind,lapply(fits,`[[`,"table")),"adjusted_HR")
   ph <- lapply(seq_along(fits),function(i){z<-capture_fit(ph_test(fits[[i]]$fit));z$exposure<-cfg$exposures[i];z})
   save_table(do.call(rbind,ph),"proportional_hazards")
 }
 if(selected("subgroups") && length(cfg$subgroups)) {
   rows <- list(); interactions <- list(); failures <- list()
   for(g in cfg$subgroups)for(v in setdiff(cfg$exposures,g)) {
     if(!is.factor(d[[g]]))stop("Subgroups must be configured factors: ",g)
     adj <- setdiff(cfg$adjustment,c(v,g))
     for(lev in levels(d[[g]])) {
       x <- d[!is.na(d[[g]]) & d[[g]]==lev,,drop=FALSE]
       ans <- tryCatch(fit_one(v,adj,x,"subgroup"),error=function(e)e)
       key <- paste(g,v,lev,sep="/")
       if(inherits(ans,"error")) failures[[key]]<-data.frame(subgroup=g,exposure=v,level=lev,error=conditionMessage(ans))
       else {z<-ans$table;z$subgroup<-g;z$level<-lev;rows[[key]]<-z}
     }
     cols <- unique(c("id","time","status",v,g,adj))
     x <- d[complete.cases(d[,cols,drop=FALSE]),cols,drop=FALSE]
     ans <- tryCatch({
       reduced<-capture_fit(cox_model(x,paste(quote_name(c(v,g,adj)),collapse=" + ")))
       rhs<-paste(c(paste(quote_name(v),quote_name(g),sep=" * "),quote_name(adj)),collapse=" + ")
       full<-capture_fit(cox_model(x,rhs));z<-compare_nested(full,reduced);z$subgroup<-g;z$exposure<-v;z
     },error=function(e)data.frame(subgroup=g,exposure=v,error=conditionMessage(e)))
     if("error" %in% names(ans))failures[[paste(g,v,"interaction")]]<-data.frame(subgroup=g,exposure=v,level="interaction",error=ans$error)
     else interactions[[paste(g,v)]]<-ans
   }
   if(length(rows))save_table(do.call(rbind,rows),"subgroup_HR")
   if(length(interactions))save_table(do.call(rbind,interactions),"interaction_tests")
   save_table(if(length(failures))do.call(rbind,failures)else data.frame(subgroup=character(),exposure=character(),level=character(),error=character()),"subgroup_failures")
 }
 if(selected("screen_plot")) {
   z<-read.csv(file.path(out,"screening_HR.csv"))
   z<-z[is.finite(z$HR)&is.finite(z$lower)&is.finite(z$upper)&z$lower>0,,drop=FALSE]
   if(!nrow(z))stop("No finite positive confidence intervals for plotting")
   pdf(file.path(out,"screening_HR.pdf"),width=9,height=max(4,0.3*nrow(z)+2))
   tryCatch({par(mar=c(5,12,2,1));plot(z$HR,seq_len(nrow(z)),log="x",xlim=range(c(z$lower,z$upper,1)),yaxt="n",ylab="",xlab="Hazard ratio (95% CI)");segments(z$lower,seq_len(nrow(z)),z$upper,seq_len(nrow(z)));axis(2,at=seq_len(nrow(z)),labels=z$term,las=1);abline(v=1,lty=2)},finally=dev.off())
 }
 save_table(data.frame(warning=warnings),paste0("warnings_",stage))
 dput(cfg,file=file.path(out,"model_specification.R"))
 capture.output(sessionInfo(),file=file.path(out,"session_info.txt"))
 invisible(out)
}
