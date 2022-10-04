library(tidyr)
library(dplyr)

final_names <- c("relief_rating#0_Coef", "relief_rating#0_Tstat", "onset_csm#0_Coef", "onset_csm#0_Tstat",
           "onset_csav#0_Coef", "onset_csav#0_Tstat", "onset_csunav#0_Coef", "onset_csunav#0_Tstat",
	   "shock#0_Coef", "shock#0_Tstat", "press#0_Coef", "press#0_Tstat")

subjects <- c("sub-FG01","sub-FG02","sub-FG03","sub-FG04","sub-FG05","sub-FG06","sub-FG07","sub-FG09","sub-FG10","sub-FG11","sub-FG12","sub-FG14","sub-FG15","sub-FG16","sub-FG17","sub-FG20","sub-FG22","sub-FG23","sub-FG25","sub-FG28R","sub-FG29R","sub-FG30R","sub-FG31R","sub-FG32R","sub-FG33R","sub-RG26","sub-RG27","sub-RG30","sub-RG31","sub-RG32","sub-RG33","sub-RG34","sub-RG35","sub-RG37","sub-RG38","sub-RG41","sub-RG44","sub-RG47","sub-RG48","sub-RG49","sub-RG50","sub-RG51","sub-RG55R","sub-RG56R","sub-RG58R", "sub-RG59R","sub-RG61R","sub-RG62R","sub-RG63R","sub-RG64R")
ROIs <- c("NAcc", "VTA", "VmPFC", "VmPFC_box")

read.subject.roi <- function(subjroi) {
   subj <- subjroi[[1]]
   roi <- subjroi[[2]]

   event_file =paste0(subj,"/func/", subj, "_task-Avoidance_events.tsv")
   X <- read.table(event_file, sep="\t", header=T, na.strings="n/a")
   omission_csav <- sort(X[X$trial_type == "omission_csav", "trial_number_cs"])
   omission_csm <- sort(X[X$trial_type == "omission_csm", "trial_number_cs"])

   csm_coefs <- paste0(rep(paste0("csm#", omission_csm), each=2), c("_Coef", "_Tstat"))
   csav_coefs <- paste0(rep(paste0("csav#", omission_csav), each=2), c("_Coef", "_Tstat"))
   names <- c("Full_Fstat", csm_coefs, csav_coefs, final_names)

   filename <- paste0("derivatives/afni/",subj,"/Avoidance/FL_results_ROI_WM_non_pooled/betas_ROI_",roi,"_REML.1D")
   Y <- read.table(filename, col.names=names)
   Y %>% tidyr::pivot_longer(tidyr::everything()) %>%
   	 tidyr::extract(name, 
   	 	"([[:alpha:]_]+).?([[0-9]+)?_([[:alpha:]]+)", 
   	 	into=c("cs", "trial_number", "type"), 
   	 	convert=T) %>%
 	 dplyr::filter(!is.na(trial_number), cs %in% c("csav", "csm")) %>%
 	 dplyr::mutate(
 	 	cs=factor(cs, levels=c("csav", "csm"), labels=c("CSav", "CSM")),
 	 	SUBJ_ID=sub("^sub-", "", subj), roi=roi, phase="AVO") %>%
 	 tidyr::pivot_wider(names_from="type")
}

X <- do.call(rbind, apply(expand.grid(subjects, ROIs), 1, FUN=read.subject.roi))
X <- X %>% tidyr::pivot_wider(names_from="roi", values_from=c("Coef", "Tstat"))
write.table(X, file="derivatives/afni/SL/SL_TbyT_WM_ROI_AVO.csv", sep=",", row.names=FALSE)
