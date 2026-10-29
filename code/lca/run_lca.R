# Categorical latent class analysis of supplied indicator matrices.
suppressPackageStartupMessages({library(poLCA);library(clue)})
args <- commandArgs(trailingOnly=TRUE)
f <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(dirname(dirname(normalizePath(f))))
cohort <- if(length(args))args[1] else "clinical"
if(!cohort %in% c("clinical","ukb"))stop("Choose clinical or ukb")
data_dir <- Sys.getenv("CEBP_DATA_DIR",file.path(root,"data"))
out <- file.path(Sys.getenv("CEBP_RESULTS_DIR",file.path(root,"results")),paste0("lca_",cohort));dir.create(out,recursive=TRUE,showWarnings=FALSE)
input <- file.path(data_dir,paste0("lca_",cohort,".csv"))
x <- read.csv(input,check.names=FALSE)
if("id" %in% names(x)){ids<-x$id;x$id<-NULL}else ids<-seq_len(nrow(x))
if(anyDuplicated(ids))stop("Duplicate IDs")
expected <- if(cohort=="clinical")c("smoking","sedentary","falls","hypertension","diabetes","glucose","haemoglobin","creatinine")else c("alcohol","smoking","fev1","fvc","grip","hba1c","tyg","glucose","cystatin","crp")
if(!identical(names(x),expected))stop("Indicator names/order must match docs/lca_input.md; outcomes must not be included")
if(anyNA(x)||!all(vapply(x,function(v)all(v %in% c(0,1)),logical(1))))stop("Provide complete binary indicators; no silent imputation/recode")
if(any(vapply(x,function(v)length(unique(v))<2,logical(1))))stop("Constant indicator detected; review coding")
write.csv(data.frame(indicator=names(x),adverse_prevalence=colMeans(x)),file.path(out,"indicator_check.csv"),row.names=FALSE)
x<-x+1L
formula<-as.formula(paste0("cbind(",paste(names(x),collapse=","),")~1"))
nrep<-as.integer(Sys.getenv("LCA_NREP","20"));maxit<-as.integer(Sys.getenv("LCA_MAXITER","3000"));B<-as.integer(Sys.getenv("LCA_BOOTSTRAPS","100"));bootrep<-as.integer(Sys.getenv("LCA_BOOT_NREP","5"))
if(nrep<1||B<1||bootrep<1||maxit<1)stop("Invalid run settings")
seed<-20261029L;selected_k<-if(cohort=="clinical")3L else 4L
# Selected K follows the manuscript specification, not an assertion that BIC chooses it.
entropy<-function(p)1+sum(p*log(pmax(p,1e-15)))/(nrow(p)*log(ncol(p)))
profiles<-function(f)sapply(f$probs,function(v)v[,2])
fit_one<-function(dat,k,reps,s){set.seed(s);w<-character();fit<-withCallingHandlers(poLCA(formula,dat,nclass=k,nrep=reps,maxiter=maxit,verbose=FALSE,calc.se=FALSE),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart("muffleWarning")});list(fit=fit,warnings=paste(unique(w),collapse="; "))}
models<-list();selection<-list()
for(k in 2:6){
 if(k*(ncol(x)+1)-1 >= 2^ncol(x)-1)stop("Too many parameters for indicator response space")
 a<-fit_one(x,k,nrep,seed+k);m<-a$fit;models[[as.character(k)]]<-m
 counts<-tabulate(m$predclass,nbins=k)
 selection[[length(selection)+1]]<-data.frame(K=k,N=nrow(x),logLik=m$llik,AIC=m$aic,BIC=m$bic,entropy=entropy(m$posterior),min_class_N=min(counts),min_class_fraction=min(counts)/nrow(x),iterations=m$numiter,at_iteration_limit=m$numiter>=maxit,warning=a$warnings)
}
sel<-do.call(rbind,selection);write.csv(sel,file.path(out,"model_selection_2_to_6.csv"),row.names=FALSE)
saveRDS(models,file.path(out,"fitted_models.rds"))
m<-models[[as.character(selected_k)]]
if(m$numiter>=maxit)stop("Selected model at iteration limit; comparison saved, bootstrap not started")
p<-profiles(m);post<-m$posterior;cl<-m$predclass
write.csv(data.frame(class=rep(seq_len(selected_k),each=ncol(p)),indicator=rep(colnames(p),selected_k),probability=as.vector(t(p))),file.path(out,"class_profiles.csv"),row.names=FALSE)
maxp<-apply(post,1,max);margin<-apply(post,1,function(v){v<-sort(v,decreasing=TRUE);v[1]-v[2]})
write.csv(data.frame(id=ids,class=cl,max_posterior=maxp,posterior_margin=margin),file.path(out,"individual_assignments.csv"),row.names=FALSE)
certainty<-do.call(rbind,lapply(seq_len(selected_k),function(g){v<-maxp[cl==g];data.frame(class=g,N=length(v),mean=mean(v),median=median(v),p10=unname(quantile(v,.1)),p90=unname(quantile(v,.9)),fraction_ge_080=mean(v>=.8),mean_margin=mean(margin[cl==g]))}))
write.csv(certainty,file.path(out,"posterior_certainty.csv"),row.names=FALSE)
boot<-list();diagnostics<-list()
for(i in seq_len(B)){
 set.seed(seed+10000L+i);ii<-sample.int(nrow(x),replace=TRUE)
 a<-tryCatch(fit_one(x[ii,],selected_k,bootrep,seed+20000L+i),error=function(e)e)
 if(inherits(a,"error")){diagnostics[[i]]<-data.frame(replicate=i,ok=FALSE,reason=conditionMessage(a));next}
 q<-profiles(a$fit)
 if(a$fit$numiter>=maxit||!all(is.finite(q))){diagnostics[[i]]<-data.frame(replicate=i,ok=FALSE,reason="Iteration limit or invalid probabilities");next}
 cost<-outer(seq_len(selected_k),seq_len(selected_k),Vectorize(function(g,h)mean((p[g,]-q[h,])^2)))
 match<-as.integer(solve_LSAP(cost));q<-q[match,,drop=FALSE]
 boot[[length(boot)+1]]<-data.frame(replicate=i,class=seq_len(selected_k),RMSE=sqrt(rowMeans((p-q)^2)),correlation=vapply(seq_len(selected_k),function(g)if(sd(p[g,])>0&&sd(q[g,])>0)cor(p[g,],q[g,])else NA_real_,numeric(1)))
 diagnostics[[i]]<-data.frame(replicate=i,ok=TRUE,reason=a$warnings)
}
write.csv(do.call(rbind,diagnostics),file.path(out,"bootstrap_diagnostics.csv"),row.names=FALSE)
if(length(boot)){
 z<-do.call(rbind,boot);write.csv(z,file.path(out,"bootstrap_replicates.csv"),row.names=FALSE)
 summary<-do.call(rbind,lapply(split(z,z$class),function(z)data.frame(class=z$class[1],successful_replicates=nrow(z),mean_RMSE=mean(z$RMSE),median_correlation=median(z$correlation,na.rm=TRUE),p025_correlation=quantile(z$correlation,.025,na.rm=TRUE),p975_correlation=quantile(z$correlation,.975,na.rm=TRUE))))
 write.csv(summary,file.path(out,"bootstrap_stability.csv"),row.names=FALSE)
}
writeLines(c(paste("Selected K:",selected_k),paste("Minimum-BIC K:",sel$K[which.min(sel$BIC)]),paste("Starts:",nrep,"Bootstrap:",B,"Bootstrap starts:",bootrep),"Class numbers are arbitrary; clinical labels require profile review."),file.path(out,"run_notes.txt"))
capture.output(sessionInfo(),file=file.path(out,"session_info.txt"))
