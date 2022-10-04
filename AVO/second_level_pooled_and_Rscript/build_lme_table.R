library(tidyr)
library(dplyr)

coefs <- c(
  "csav#0_Coef",
  "csm#0_Coef",
  "onset_csm#0_Coef",
  "onset_csav#0_Coef",
  "onset_csunav#0_Coef"
)

subjects <- c("sub-FG01","sub-FG02","sub-FG03","sub-FG04","sub-FG05","sub-FG06","sub-FG07","sub-FG09","sub-FG10","sub-FG11","sub-FG12","sub-FG14","sub-FG15","sub-FG16","sub-FG17","sub-FG20","sub-FG22","sub-FG23","sub-FG25","sub-FG28R","sub-FG29R","sub-FG30R","sub-FG31R","sub-FG32R","sub-FG33R","sub-RG26","sub-RG27","sub-RG30","sub-RG31","sub-RG32","sub-RG33","sub-RG34","sub-RG35","sub-RG37","sub-RG38","sub-RG41","sub-RG44","sub-RG47","sub-RG48","sub-RG49","sub-RG50","sub-RG51","sub-RG55R","sub-RG56R","sub-RG58R", "sub-RG59R","sub-RG61R","sub-RG62R","sub-RG63R","sub-RG64R")

do.subject <- function(subj) {
   filenames <- paste(paste0("derivatives/afni/",subj,"/Avoidance/FL_results_WB_pooled/betas_REML+tlrc"), "'[", coefs, "]'", sep="")
   data.frame(Subj=subj, Phase=factor(c('OM', 'OM', 'AN', 'AN', 'AN')), CS=factor(c('csav', 'csm', 'csm', 'csav', 'csunav')), InputFile=filenames)
}

X <- do.call(rbind, lapply(subjects, FUN=do.subject))
print(X)
write.table(X, file="derivatives/afni/SL/AVO/lme_table.csv", sep="\t", row.names=FALSE, quote=FALSE)
