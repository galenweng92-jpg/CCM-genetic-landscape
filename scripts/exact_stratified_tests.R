#!/usr/bin/env Rscript
# Base-R statistical cross-check. No continuity correction for primary CMH result;
# corrected CMH is supplied separately. Zero-cell strata are not pseudo-counted.
args <- commandArgs(trailingOnly=TRUE)
out <- args[[1]]
counts <- read.csv(file.path(out,"association_2x2_counts.csv"),fileEncoding="UTF-8-BOM",check.names=FALSE)
strata <- read.csv(file.path(out,"association_study_strata.csv"),fileEncoding="UTF-8-BOM",check.names=FALSE)
results <- lapply(seq_len(nrow(counts)),function(i){
  r <- counts[i,]
  t <- matrix(c(r$a_driver_pos_ccm_pos,r$b_driver_pos_ccm_neg,
                r$c_driver_neg_ccm_pos,r$d_driver_neg_ccm_neg),nrow=2,byrow=TRUE)
  f <- fisher.test(t,alternative="two.sided",conf.level=.95)
  ss <- strata[strata$analysis_id==r$analysis_id,]
  s <- ss[tolower(as.character(ss$informative_for_CMH))=="true",]
  arr <- array(NA_real_,dim=c(2,2,nrow(s)))
  for(j in seq_len(nrow(s))) arr[,,j] <- matrix(unlist(s[j,c("a","b","c","d")]),2,byrow=TRUE)
  mh <- mantelhaen.test(arr,correct=FALSE,conf.level=.95)
  mhc <- mantelhaen.test(arr,correct=TRUE,conf.level=.95)
  mhe <- mantelhaen.test(arr,exact=TRUE,conf.level=.95)
  data.frame(r,
             odds_ratio_crude_sample=if(t[1,2]*t[2,1]>0)t[1,1]*t[2,2]/(t[1,2]*t[2,1])else NA_real_,
             fisher_conditional_odds_ratio=unname(f$estimate),
             fisher_exact_ci95_lower=f$conf.int[1],
             fisher_exact_ci95_upper=f$conf.int[2],
             fisher_two_sided_p=f$p.value,
             strata_total=nrow(ss),
             strata_informative=nrow(s),
             informative_n=sum(s$n),
             mh_odds_ratio=unname(mh$estimate),
             mh_ci95_lower=mh$conf.int[1],
             mh_ci95_upper=mh$conf.int[2],
             cmh_chisq=unname(mh$statistic),
             cmh_p_uncorrected=mh$p.value,
             cmh_p_continuity_corrected=mhc$p.value,
             stratified_exact_conditional_or=unname(mhe$estimate),
             stratified_exact_ci95_lower=mhe$conf.int[1],
             stratified_exact_ci95_upper=mhe$conf.int[2],
             stratified_exact_two_sided_p=mhe$p.value,
             mh_asymptotic_ci_estimable=all(is.finite(mh$conf.int)))
})
results <- do.call(rbind,results)
write.csv(results,file.path(out,"association_exact_stratified_results.csv"),row.names=FALSE,na="")
jsonlite::write_json(list(R_version=R.version.string,methods=list(
  Fisher="stats::fisher.test, two-sided, conditional maximum-likelihood OR with exact CI",
  MH="stats::mantelhaen.test; common Mantel-Haenszel OR with asymptotic log-scale 95% CI",
  CMH="Two-sided asymptotic chi-square test, primary without continuity correction; corrected p also supplied",
  exact_stratified="stats::mantelhaen.test(exact=TRUE); conditional exact common OR, CI and two-sided p; use for sparse all-three-MAP3K3 comparison where asymptotic MH CI is undefined",
  strata="Study / report field; strata with no variation in exposure or outcome excluded as non-informative",
  independence="Exploratory source-record inference only, not proof of molecular exclusivity or patient-level risk"
),results=results),file.path(out,"association_exact_stratified_results.json"),pretty=TRUE,digits=16)
print(results[,c("analysis_id","mh_odds_ratio","mh_ci95_lower","mh_ci95_upper","fisher_two_sided_p","cmh_p_uncorrected")],row.names=FALSE)

# Independent Wilson cross-check against stats::prop.test(correct=FALSE).
rate_files <- c("familial_composition_wilson.csv","paired_endpoints_wilson.csv",
                "sporadic_detection_sensitivity_wilson.csv",
                "map3k3_pik3ca_cofinding_sensitivity_wilson.csv",
                "study_specific_detection_wilson.csv")
checked <- 0
for (fn in rate_files) {
  rr <- read.csv(file.path(out,fn),fileEncoding="UTF-8-BOM")
  for (i in seq_len(nrow(rr))) {
    if (rr$eligible_n[i]==0) next
    ci <- suppressWarnings(prop.test(rr$positive[i],rr$eligible_n[i],correct=FALSE)$conf.int)
    stopifnot(max(abs(ci-c(rr$ci95_lower[i],rr$ci95_upper[i])))<1e-10)
    checked <- checked+1
  }
}
jsonlite::write_json(list(wilson_rows_checked=checked,reference="stats::prop.test(correct=FALSE)",
                          maximum_tolerance=1e-10,passed=TRUE),
                     file.path(out,"wilson_independent_R_validation.json"),pretty=TRUE)
