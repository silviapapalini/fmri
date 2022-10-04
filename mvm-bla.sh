#!/bin/sh

dirA=derivatives/afni
dirP=derivatives/fmriprep
results_dir=$dirA/SL/AVO
mask=$results_dir/mask_union_task-Avoidance.nii.gz
if [ ! -e $mask_dset ]; then
	3dmask_tool -input $dirP/sub-*/func/sub-*_task-Avoidance_space-MNI152NLin20011cAsym_desc-brain_mask.nii.gz -prefix $mask -union
fi

#mask=derivatives/afni/SL/mean_gm_mask_resample.nii.gz

3dMVM -prefix $results_dir/phase_cs -jobs 24 \
  -mask $mask \
  -bsVars 'Group' \
  -wsVars 'CS*Phase' \
  -dataTable                                                                   \
Subj		Group	Phase	CS	InputFile \
sub-FG01	FG	OM	csav	derivatives/afni/sub-FG01/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG01	FG	OM	csm	derivatives/afni/sub-FG01/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG01	FG	AN	csm	derivatives/afni/sub-FG01/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG01	FG	AN	csav	derivatives/afni/sub-FG01/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG01	FG	AN	csunav	derivatives/afni/sub-FG01/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG02	FG	OM	csav	derivatives/afni/sub-FG02/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG02	FG	OM	csm	derivatives/afni/sub-FG02/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG02	FG	AN	csm	derivatives/afni/sub-FG02/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG02	FG	AN	csav	derivatives/afni/sub-FG02/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG02	FG	AN	csunav	derivatives/afni/sub-FG02/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG03	FG	OM	csav	derivatives/afni/sub-FG03/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG03	FG	OM	csm	derivatives/afni/sub-FG03/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG03	FG	AN	csm	derivatives/afni/sub-FG03/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG03	FG	AN	csav	derivatives/afni/sub-FG03/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG03	FG	AN	csunav	derivatives/afni/sub-FG03/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG04	FG	OM	csav	derivatives/afni/sub-FG04/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG04	FG	OM	csm	derivatives/afni/sub-FG04/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG04	FG	AN	csm	derivatives/afni/sub-FG04/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG04	FG	AN	csav	derivatives/afni/sub-FG04/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG04	FG	AN	csunav	derivatives/afni/sub-FG04/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG05	FG	OM	csav	derivatives/afni/sub-FG05/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG05	FG	OM	csm	derivatives/afni/sub-FG05/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG05	FG	AN	csm	derivatives/afni/sub-FG05/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG05	FG	AN	csav	derivatives/afni/sub-FG05/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG05	FG	AN	csunav	derivatives/afni/sub-FG05/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG06	FG	OM	csav	derivatives/afni/sub-FG06/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG06	FG	OM	csm	derivatives/afni/sub-FG06/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG06	FG	AN	csm	derivatives/afni/sub-FG06/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG06	FG	AN	csav	derivatives/afni/sub-FG06/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG06	FG	AN	csunav	derivatives/afni/sub-FG06/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG55R	RG	OM	csav	derivatives/afni/sub-RG55R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG55R	RG	OM	csm	derivatives/afni/sub-RG55R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG55R	RG	AN	csm	derivatives/afni/sub-RG55R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG55R	RG	AN	csav	derivatives/afni/sub-RG55R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG55R	RG	AN	csunav	derivatives/afni/sub-RG55R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG56R	RG	OM	csav	derivatives/afni/sub-RG56R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG56R	RG	OM	csm	derivatives/afni/sub-RG56R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG56R	RG	AN	csm	derivatives/afni/sub-RG56R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG56R	RG	AN	csav	derivatives/afni/sub-RG56R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG56R	RG	AN	csunav	derivatives/afni/sub-RG56R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG58R	RG	OM	csav	derivatives/afni/sub-RG58R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG58R	RG	OM	csm	derivatives/afni/sub-RG58R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG58R	RG	AN	csm	derivatives/afni/sub-RG58R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG58R	RG	AN	csav	derivatives/afni/sub-RG58R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG58R	RG	AN	csunav	derivatives/afni/sub-RG58R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG59R	RG	OM	csav	derivatives/afni/sub-RG59R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG59R	RG	OM	csm	derivatives/afni/sub-RG59R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG59R	RG	AN	csm	derivatives/afni/sub-RG59R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG59R	RG	AN	csav	derivatives/afni/sub-RG59R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG59R	RG	AN	csunav	derivatives/afni/sub-RG59R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG61R	RG	OM	csav	derivatives/afni/sub-RG61R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG61R	RG	OM	csm	derivatives/afni/sub-RG61R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG61R	RG	AN	csm	derivatives/afni/sub-RG61R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG61R	RG	AN	csav	derivatives/afni/sub-RG61R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG61R	RG	AN	csunav	derivatives/afni/sub-RG61R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG62R	RG	OM	csav	derivatives/afni/sub-RG62R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG62R	RG	OM	csm	derivatives/afni/sub-RG62R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG62R	RG	AN	csm	derivatives/afni/sub-RG62R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG62R	RG	AN	csav	derivatives/afni/sub-RG62R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG62R	RG	AN	csunav	derivatives/afni/sub-RG62R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG63R	RG	OM	csav	derivatives/afni/sub-RG63R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG63R	RG	OM	csm	derivatives/afni/sub-RG63R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG63R	RG	AN	csm	derivatives/afni/sub-RG63R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG63R	RG	AN	csav	derivatives/afni/sub-RG63R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG63R	RG	AN	csunav	derivatives/afni/sub-RG63R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG64R	RG	OM	csav	derivatives/afni/sub-RG64R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG64R	RG	OM	csm	derivatives/afni/sub-RG64R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG64R	RG	AN	csm	derivatives/afni/sub-RG64R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG64R	RG	AN	csav	derivatives/afni/sub-RG64R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG64R	RG	AN	csunav	derivatives/afni/sub-RG64R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]'
