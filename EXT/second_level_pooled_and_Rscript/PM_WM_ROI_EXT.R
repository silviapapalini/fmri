library(tidyr)
library(dplyr)

final_names <- c("relief_rating#0_Coef", "relief_rating#0_Tstat", "onset_csm#0_Coef", "onset_csm#0_Tstat",
           "onset_csav#0_Coef", "onset_csav#0_Tstat", "csavnonMcsavM#0_Coef", "csavnonMcsavM#0_Tstat",
           "csmnonMcsmM#0_Coef", "csmnonMcsmM#0_Tstat",
           "csavnonMANDcsmnonMcsavMANDcsmM#0_Coef", "csavnonMANDcsmnonMcsavMANDcsmM#0_Tstat",
           "csavnonMcsmnonM#0_Coef", "csavnonMcsmnonM#0_Tstat",
           "csmMcsavM#0_Coef", "csmMcsavM#0_Tstat", "csavmodulatedcsm_modulated#0_Coef","csavmodulatedcsm_modulated#0_Tstat",
           "csav_modulatedANDcsm_modulated#0_Coef","csav_modulatedANDcsm_modulated#0_Tstat")

#0=csav+csm
#1=csav only
#2=csm only
#3=nothing

names0 <- c("csav#0_Coef", "csav#0_Tstat", "csav#1_Coef", "csav#1_Tstat", "csm#0_Coef", "csm#0_Tstat", "csm#1_Coef", "csm#1_Tstat")
names1 <- c("csav#0_Coef", "csav#0_Tstat", "csav#1_Coef", "csav#1_Tstat", "csm#0_Coef", "csm#0_Tstat")
names2 <- c("csav#0_Coef", "csav#0_Tstat", "csm#0_Coef", "csm#0_Tstat", "csm#1_Coef", "csm#1_Tstat")
names3 <- c("csav#0_Coef", "csav#0_Tstat", "csm#0_Coef", "csm#0_Tstat")

subjects_coefs <- c("sub-FG01"=0, "sub-FG02"=0, "sub-FG03"=0, "sub-FG04"=0, "sub-FG05"=3, "sub-FG06"=0, "sub-FG07"=1, "sub-FG09"=0, "sub-FG10"=0, "sub-FG11"=0, "sub-FG12"=0, "sub-FG14"=0, "sub-FG16"=0, "sub-FG17"=0, "sub-FG20"=0, "sub-FG22"=0, "sub-FG23"=0, "sub-FG25"=0, "sub-FG28R"=0, "sub-FG29R"=0, "sub-FG30R"=0, "sub-FG31R"=0, "sub-FG32R"=0, "sub-FG33R"=0, "sub-RG26"=0, "sub-RG27"=0, "sub-RG30"=1, "sub-RG31"=0, "sub-RG34"=0, "sub-RG37"=0, "sub-RG38"=0, "sub-RG41"=0, "sub-RG44"=0, "sub-RG47"=0, "sub-RG48"=0, "sub-RG49"=0, "sub-RG50"=0, "sub-RG51"=0, "sub-RG55R"=1, "sub-RG56R"=0, "sub-RG58R"=0, "sub-RG59R"=0, "sub-RG61R"=0, "sub-RG62R"=1, "sub-RG63R"=0, "sub-RG64R"=0)
ROIs <- c("VTA","canlab_VTA_bilateral", "NAcc", "canlab_NAC_bilateral","VmPFC", "VmPFC_box", "VmPFC_csp_vs_csm")

read.subject.roi <- function(subjroi) {
   subj <- subjroi[[1]]
   coefs <- subjects_coefs[subj]
   roi <- subjroi[[2]]
   tryCatch({
	   filename <- paste0("derivatives/afni/",subj,"/Extinction/FL_results_ROI_WM_parametric_modulation/betas_ROI_",roi,"_REML.1D")
	   coef_names <- switch(coefs + 1, names0, names1, names2, names3)
	   names <- c("Full_Fstat", coef_names, final_names)
	   Y <- read.table(filename, col.names=names)
	   Y %>% tidyr::pivot_longer(tidyr::everything()) %>%
	   	 tidyr::extract(name, 
	   	 	"([[:alpha:]]+).?([[0-9]+)?_([[:alpha:]]+)", 
	   	 	into=c("event_type", "count", "type"), 
	   	 	convert=T) %>%
	   	 dplyr::filter(!is.na(count)) %>%
	 	 dplyr::mutate(subj_id=sub("^sub-", "", subj), roi=roi, phase="EXT") %>%
	 	 tidyr::pivot_wider(names_from="type")
   }, error=function(e) { cat("error", subj, roi); stop(e) })
}

X <- do.call(rbind, apply(expand.grid(names(subjects_coefs), ROIs), 1, FUN=read.subject.roi))
X <- X %>% tidyr::pivot_wider(names_from="roi", values_from=c("Coef", "Tstat"))
write.table(X, file="derivatives/afni/SL/EXT/SL_PM_WM_ROI_EXT.csv", sep=",", row.names=FALSE)
