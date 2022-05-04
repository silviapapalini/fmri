library(ggplot2)
library(dplyr)

setwd("~/FASTING_MRI_BIDS/derivatives/afni/SL/")
X <- read.csv("SL_pooled_roi_EXT.csv", header=T) %>%
  tidyr::extract(SUBJ_ID, c("GROUP", "ID"), "([A-Z]+)([0-9]+R?)") %>%
  tidyr::pivot_longer(c(CSAV,CSM,CSAV_CSM,CSM_CSAV), names_to="CS") %>%
  dplyr::filter(CS=="CSM_CSAV") %>% 
  dplyr::filter(ROI=="VmPFC")
 

ggplot(X) + geom_boxplot(aes(x=GROUP, y=value, group=GROUP)) + facet_wrap(ROI~CS)
ggplot(X) + geom_jitter(aes(x=GROUP, y=value, group=GROUP, color=GROUP), width=0.1) + facet_wrap(ROI~CS)

p.values <- X %>%
  dplyr::group_by(ROI, CS, GROUP) %>%
  summarise(value = list(value)) %>%
  tidyr::pivot_wider(names_from="GROUP", values_from=value) %>%
  dplyr::summarise(p_value = t.test(unlist(FG), unlist(RG))$p.value)

