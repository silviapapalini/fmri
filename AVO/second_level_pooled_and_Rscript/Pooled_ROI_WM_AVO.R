library(tidyr)
library(dplyr)

names <- c(
  "Full_Fstat",
  "csav#0_Coef",
  "csav#0_Tstat",
  "csm#0_Coef",
  "csm#0_Tstat",
  "relief_rating#0_Coef",
  "relief_rating#0_Tstat",
  "onset_csm#0_Coef",
  "onset_csm#0_Tstat",
  "onset_csav#0_Coef",
  "onset_csav#0_Tstat",
   "onset_csunav#0_Coef",
   "onset_csunav#0_Tstat",
   "shock#0_Coef",
   "shock#0_Tstat",
   "press#0_Coef",
   "press#0_Tstat",
   "csav_vs_csm#0_Coef",
   "csav_vs_csm#0_Tstat",
   "csav_and_csm#0_Coef",
   "csav_and_csm#0_Tstat",
   "csm_vs_csav#0_Coef",
   "csm_vs_csav#0_Tstat",
   "onset_csav_vs_onset_csm#0_Coef",
   "onset_csav_vs_onset_csm#0_Tstat",
   "onset_csm_vs_onset_csav#0_Coef",
   "onset_csm_vs_onset_csav#0_Tstat",
   "onset_csm_vs_onset_csunav#0_Coef",
   "onset_csm_vs_onset_csunav#0_Tstat",
   "onset_csav_vs_onset_csunav#0_Coef",
   "onset_csav_vs_onset_csunav#0_Tstat"
)

subjects <- c("sub-FG01","sub-FG02","sub-FG03","sub-FG04","sub-FG05","sub-FG06","sub-FG07","sub-FG09","sub-FG10","sub-FG11","sub-FG12","sub-FG14","sub-FG15","sub-FG16","sub-FG17","sub-FG20","sub-FG22","sub-FG23","sub-FG25","sub-FG28R","sub-FG29R","sub-FG30R","sub-FG31R","sub-FG32R","sub-FG33R","sub-RG26","sub-RG27","sub-RG30","sub-RG31","sub-RG32","sub-RG33","sub-RG34","sub-RG35","sub-RG37","sub-RG38","sub-RG41","sub-RG44","sub-RG47","sub-RG48","sub-RG49","sub-RG50","sub-RG51","sub-RG55R","sub-RG56R","sub-RG58R", "sub-RG59R","sub-RG61R","sub-RG62R","sub-RG63R","sub-RG64R")
ROIs <- c("VTA","canlab_VTA_bilateral", "NAcc", "canlab_NAC_bilateral","VmPFC", "VmPFC_box", "VmPFC_csp_vs_csm")

read.subject.roi <- function(subjroi) {
   subj <- subjroi[[1]]
   roi <- subjroi[[2]]
   filename <- paste0("derivatives/afni/",subj,"/Avoidance/FL_results_ROI_WM_pooled/betas_ROI_",roi,"_REML.1D")
   Y <- read.table(filename, col.names=names)
   Y %>% tidyr::pivot_longer(tidyr::everything()) %>%
   	 tidyr::extract(name, 
   	 	"([[:alpha:]_]+).?([[0-9]+)?_([[:alpha:]]+)$", 
   	 	into=c("event_type", "count", "type"), 
   	 	convert=T) %>%
 	 dplyr::filter(!is.na(count)) %>%
 	 dplyr::mutate(SUBJ_ID=sub("^sub-", "", subj), roi=roi, phase="AVO") %>%
 	 tidyr::pivot_wider(names_from="type")
}

X <- do.call(rbind, apply(expand.grid(subjects, ROIs), 1, FUN=read.subject.roi))
print(X)
X <- X %>% tidyr::pivot_wider(names_from="roi", values_from=c("Coef", "Tstat"))
write.table(X, file="derivatives/afni/SL/AVO/SL_pooled_ROI_WM_AVO.csv", sep=",", row.names=FALSE)
