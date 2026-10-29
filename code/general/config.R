# Edit variables and time units to match the approved analysis specification.
config <- list(input="analysis.csv",factors=c("phenotype","age_category","education","adt","chemotherapy"),
 indicators=character(),cox_rhs="score + age_category + education + adt + chemotherapy",
 fg_rhs="score + age_category + education + adt + chemotherapy",
 screen=character(),screen_adjust=character(),group="phenotype",time_unit="years",
 period_cuts=c(3,6),lags=c(0,1,2),model_rhs=list(base="age_category + education + adt + chemotherapy",score="age_category + education + adt + chemotherapy + score",phenotype="age_category + education + adt + chemotherapy + phenotype",combined="age_category + education + adt + chemotherapy + score + phenotype"))
