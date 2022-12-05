library(tidyr)
library(dplyr)

names <- c(
  "Full_Fstat",
  "csav_early#0_Coef",
  "csav_early#0_Tstat",
  "csm_early#0_Coef",
  "csm_early#0_Tstat",
  "csav_late#0_Coef",
  "csav_late#0_Tstat",
  "csm_late#0_Coef",
  "csm_late#0_Tstat",
  "relief_rating#0_Coef",
  "relief_rating#0_Tstat",
  "onset_csm_early#0_Coef",
  "onset_csm_early#0_Tstat",
  "onset_csav_early#0_Coef",
  "onset_csav_early#0_Tstat",
  "onset_csm_late#0_Coef",
  "onset_csm_late#0_Tstat",
  "onset_csav_late#0_Coef",
  "onset_csav_late#0_Tstat", 
   "csav_early_vs_csm_early#0_Coef",
   "csav_early_vs_csm_early#0_Tstat",
   "csm_early_vs_csav_early#0_Coef",
   "csm_early_vs_csav_early#0_Tstat",
   "csav_late_vs_csm_late#0_Coef",
   "csav_late_vs_csm_late#0_Tstat",
   "csm_late_vs_csav_late#0_Coef",
   "csm_late_vs_csav_late#0_Tstat",
   "onset_csav_early_vs_onset_csm_early#0_Coef",
   "onset_csav_early_vs_onset_csm_early#0_Tstat",
   "onset_csm_early_vs_onset_csav_early#0_Coef",
   "onset_csm_early_vs_onset_csav_early#0_Tstat",
   "onset_csav_late_vs_onset_csm_late#0_Coef",
   "onset_csav_late_vs_onset_csm_late#0_Tstat",
   "onset_csm_late_vs_onset_csav_late#0_Coef",
   "onset_csm_late_vs_onset_csav_late#0_Tstat",
   "onset_csav_early_vs_late_vs_onset_csm_early_vs_onset_csm_late#0_Coef", 
   "onset_csav_early_vs_late_vs_onset_csm_early_vs_onset_csm_late#0_Tstat", 	    
   "csav_early_AND_csm_early_vs_csav_late_AND_csm_late#0_Coef",
   "csav_early_AND_csm_early_vs_csav_late_AND_csm_late#0_Tstat", 	    
   "onset_csav_early_AND_onset_csm_early_vs_onset_csav_late_AND_onset_csm_late#0_Coef", 
   "onset_csav_early_AND_onset_csm_early_vs_onset_csav_late_AND_onset_csm_late#0_Tstat", 	   
   "csm_early_vs_csav_early_minus_csm_late_vs_csav_late#0_Coef",
   "csm_early_vs_csav_early_minus_csm_late_vs_csav_late#0_Tstat",
   "csav_early_vs_csm_early_minus_csav_late_vs_csm_late#0_Coef", 
   "csav_early_vs_csm_early_minus_csav_late_vs_csm_late#0__Tstat" 
)

subjects <- c("sub-FG01","sub-FG02","sub-FG03","sub-FG04","sub-FG05","sub-FG06","sub-FG07","sub-FG09","sub-FG10","sub-FG11","sub-FG12","sub-FG14","sub-FG16","sub-FG17","sub-FG20","sub-FG22","sub-FG23","sub-FG25","sub-FG28R","sub-FG29R","sub-FG30R","sub-FG31R","sub-FG32R","sub-FG33R","sub-RG26","sub-RG27","sub-RG30","sub-RG31","sub-RG34","sub-RG37","sub-RG38","sub-RG41","sub-RG44","sub-RG47","sub-RG48","sub-RG49","sub-RG50","sub-RG51","sub-RG55R","sub-RG56R","sub-RG58R","sub-RG59R","sub-RG61R","sub-RG62R","sub-RG63R","sub-RG64R")

ROIs <- c("VTA","canlab_VTA_bilateral", "NAcc", "canlab_NAC_bilateral","VmPFC", "VmPFC_box", "VmPFC_csp_vs_csm")

read.subject.roi <- function(subjroi) {
   subj <- subjroi[[1]]
   roi <- subjroi[[2]]
   filename <- paste0("derivatives/afni/",subj,"/Extinction/FL_results_ROI_WM_pooled_early_late/betas_ROI_",roi,"_REML.1D")
   Y <- read.table(filename, col.names=names)
   Y %>% tidyr::pivot_longer(tidyr::everything()) %>%
   	 tidyr::extract(name, 
   	 	"([[:alpha:]_]+).?([[0-9]+)?_([[:alpha:]]+)$", 
   	 	into=c("event_type", "count", "type"), 
   	 	convert=T) %>%
 	 dplyr::filter(!is.na(count)) %>%
 	 dplyr::mutate(SUBJ_ID=sub("^sub-", "", subj), roi=roi, phase="EXT") %>%
 	 tidyr::pivot_wider(names_from="type")
}

X <- do.call(rbind, apply(expand.grid(subjects, ROIs), 1, FUN=read.subject.roi))
print(X)
X <- X %>% tidyr::pivot_wider(names_from="roi", values_from=c("Coef", "Tstat"))
write.table(X, file="derivatives/afni/SL/EXT/SL_pooled_ROI_WM_EXT_early_late.csv", sep=",", row.names=FALSE)
