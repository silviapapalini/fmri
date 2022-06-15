library(tidyr)
library(dplyr)

#X <- read.table("sub-FG01/func/sub-FG01_task-Extinction_events.tsv", header=T, sep="\t")
#X <- X[with(X, order(trial_type, onset)),]
#X[X$trial_type %in% c("relief csav", "relief csm"), "trial_number"] <- rep(1:8, times=2)

names <- readLines("~/fmri-process/EXT/sub-bricks-from-FL_EXT")

subjects <- c("sub-FG01","sub-FG02","sub-FG03","sub-FG04","sub-FG05","sub-FG06","sub-FG07","sub-FG09","sub-FG10","sub-FG11","sub-FG12","sub-FG14","sub-FG16","sub-FG17","sub-FG20","sub-FG22","sub-FG23","sub-FG25","sub-FG28R","sub-FG29R","sub-FG30R","sub-FG31R","sub-FG32R","sub-FG33R","sub-RG26","sub-RG27","sub-RG30","sub-RG31","sub-RG34","sub-RG37","sub-RG38","sub-RG41","sub-RG44","sub-RG47","sub-RG48","sub-RG49","sub-RG50","sub-RG51","sub-RG55R","sub-RG56R","sub-RG58R","sub-RG59R","sub-RG61R","sub-RG62R","sub-RG63R","sub-RG64R")
ROIs <- c("NAcc", "VTA", "VmPFC")

read.subject.roi <- function(subjroi) {
   subj <- subjroi[[1]]
   roi <- subjroi[[2]]
   filename <- paste0("derivatives/afni/",subj,"/Extinction/FL_results_ROI_WM_non_pooled/betas_ROI_",roi,"_REML.1D")
   Y <- read.table(filename, col.names=names)
   Y %>% tidyr::pivot_longer(tidyr::everything()) %>%
   	 tidyr::extract(name, 
   	 	"([[:alpha:]]+).?([[0-9])?_([[:alpha:]]+)", 
   	 	into=c("CS", "trial_number", "type"), 
   	 	convert=T) %>%
 	 dplyr::filter(!is.na(trial_number)) %>%
 	 dplyr::mutate(
 	 	CS=factor(CS, levels=c("csav", "csm"), labels=c("CSav", "CSM")),
 	 	SUBJ_ID=sub("^sub-", "", subj), roi=roi, X.TRIAL_POOL="_3", phase="EXT") %>%
 	 tidyr::pivot_wider(names_from="type")
}

X <- do.call(rbind, apply(expand.grid(subjects, ROIs), 1, FUN=read.subject.roi))
X <- X %>% tidyr::pivot_wider(names_from="roi", values_from=c("Coef", "Tstat"))
write.table(X, file="derivatives/afni/SL/EXT/SL_TbyT_ROI_WM_EXT.csv", sep=",", row.names=FALSE)
