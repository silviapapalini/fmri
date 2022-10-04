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
sub-FG07	FG	OM	csav	derivatives/afni/sub-FG07/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG07	FG	OM	csm	derivatives/afni/sub-FG07/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG07	FG	AN	csm	derivatives/afni/sub-FG07/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG07	FG	AN	csav	derivatives/afni/sub-FG07/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG07	FG	AN	csunav	derivatives/afni/sub-FG07/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG09	FG	OM	csav	derivatives/afni/sub-FG09/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG09	FG	OM	csm	derivatives/afni/sub-FG09/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG09	FG	AN	csm	derivatives/afni/sub-FG09/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG09	FG	AN	csav	derivatives/afni/sub-FG09/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG09	FG	AN	csunav	derivatives/afni/sub-FG09/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG10	FG	OM	csav	derivatives/afni/sub-FG10/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG10	FG	OM	csm	derivatives/afni/sub-FG10/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG10	FG	AN	csm	derivatives/afni/sub-FG10/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG10	FG	AN	csav	derivatives/afni/sub-FG10/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG10	FG	AN	csunav	derivatives/afni/sub-FG10/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG11	FG	OM	csav	derivatives/afni/sub-FG11/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG11	FG	OM	csm	derivatives/afni/sub-FG11/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG11	FG	AN	csm	derivatives/afni/sub-FG11/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG11	FG	AN	csav	derivatives/afni/sub-FG11/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG11	FG	AN	csunav	derivatives/afni/sub-FG11/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG12	FG	OM	csav	derivatives/afni/sub-FG12/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG12	FG	OM	csm	derivatives/afni/sub-FG12/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG12	FG	AN	csm	derivatives/afni/sub-FG12/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG12	FG	AN	csav	derivatives/afni/sub-FG12/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG12	FG	AN	csunav	derivatives/afni/sub-FG12/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG14	FG	OM	csav	derivatives/afni/sub-FG14/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG14	FG	OM	csm	derivatives/afni/sub-FG14/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG14	FG	AN	csm	derivatives/afni/sub-FG14/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG14	FG	AN	csav	derivatives/afni/sub-FG14/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG14	FG	AN	csunav	derivatives/afni/sub-FG14/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG15	FG	OM	csav	derivatives/afni/sub-FG15/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG15	FG	OM	csm	derivatives/afni/sub-FG15/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG15	FG	AN	csm	derivatives/afni/sub-FG15/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG15	FG	AN	csav	derivatives/afni/sub-FG15/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG15	FG	AN	csunav	derivatives/afni/sub-FG15/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG16	FG	OM	csav	derivatives/afni/sub-FG16/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG16	FG	OM	csm	derivatives/afni/sub-FG16/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG16	FG	AN	csm	derivatives/afni/sub-FG16/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG16	FG	AN	csav	derivatives/afni/sub-FG16/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG16	FG	AN	csunav	derivatives/afni/sub-FG16/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG17	FG	OM	csav	derivatives/afni/sub-FG17/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG17	FG	OM	csm	derivatives/afni/sub-FG17/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG17	FG	AN	csm	derivatives/afni/sub-FG17/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG17	FG	AN	csav	derivatives/afni/sub-FG17/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG17	FG	AN	csunav	derivatives/afni/sub-FG17/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG20	FG	OM	csav	derivatives/afni/sub-FG20/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG20	FG	OM	csm	derivatives/afni/sub-FG20/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG20	FG	AN	csm	derivatives/afni/sub-FG20/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG20	FG	AN	csav	derivatives/afni/sub-FG20/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG20	FG	AN	csunav	derivatives/afni/sub-FG20/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG22	FG	OM	csav	derivatives/afni/sub-FG22/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG22	FG	OM	csm	derivatives/afni/sub-FG22/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG22	FG	AN	csm	derivatives/afni/sub-FG22/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG22	FG	AN	csav	derivatives/afni/sub-FG22/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG22	FG	AN	csunav	derivatives/afni/sub-FG22/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG23	FG	OM	csav	derivatives/afni/sub-FG23/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG23	FG	OM	csm	derivatives/afni/sub-FG23/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG23	FG	AN	csm	derivatives/afni/sub-FG23/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG23	FG	AN	csav	derivatives/afni/sub-FG23/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG23	FG	AN	csunav	derivatives/afni/sub-FG23/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG25	FG	OM	csav	derivatives/afni/sub-FG25/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG25	FG	OM	csm	derivatives/afni/sub-FG25/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG25	FG	AN	csm	derivatives/afni/sub-FG25/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG25	FG	AN	csav	derivatives/afni/sub-FG25/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG25	FG	AN	csunav	derivatives/afni/sub-FG25/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG28R	FG	OM	csav	derivatives/afni/sub-FG28R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG28R	FG	OM	csm	derivatives/afni/sub-FG28R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG28R	FG	AN	csm	derivatives/afni/sub-FG28R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG28R	FG	AN	csav	derivatives/afni/sub-FG28R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG28R	FG	AN	csunav	derivatives/afni/sub-FG28R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG29R	FG	OM	csav	derivatives/afni/sub-FG29R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG29R	FG	OM	csm	derivatives/afni/sub-FG29R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG29R	FG	AN	csm	derivatives/afni/sub-FG29R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG29R	FG	AN	csav	derivatives/afni/sub-FG29R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG29R	FG	AN	csunav	derivatives/afni/sub-FG29R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG30R	FG	OM	csav	derivatives/afni/sub-FG30R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG30R	FG	OM	csm	derivatives/afni/sub-FG30R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG30R	FG	AN	csm	derivatives/afni/sub-FG30R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG30R	FG	AN	csav	derivatives/afni/sub-FG30R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG30R	FG	AN	csunav	derivatives/afni/sub-FG30R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG31R	FG	OM	csav	derivatives/afni/sub-FG31R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG31R	FG	OM	csm	derivatives/afni/sub-FG31R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG31R	FG	AN	csm	derivatives/afni/sub-FG31R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG31R	FG	AN	csav	derivatives/afni/sub-FG31R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG31R	FG	AN	csunav	derivatives/afni/sub-FG31R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG32R	FG	OM	csav	derivatives/afni/sub-FG32R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG32R	FG	OM	csm	derivatives/afni/sub-FG32R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG32R	FG	AN	csm	derivatives/afni/sub-FG32R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG32R	FG	AN	csav	derivatives/afni/sub-FG32R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG32R	FG	AN	csunav	derivatives/afni/sub-FG32R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-FG33R	FG	OM	csav	derivatives/afni/sub-FG33R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-FG33R	FG	OM	csm	derivatives/afni/sub-FG33R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-FG33R	FG	AN	csm	derivatives/afni/sub-FG33R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-FG33R	FG	AN	csav	derivatives/afni/sub-FG33R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-FG33R	FG	AN	csunav	derivatives/afni/sub-FG33R/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG26	RG	OM	csav	derivatives/afni/sub-RG26/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG26	RG	OM	csm	derivatives/afni/sub-RG26/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG26	RG	AN	csm	derivatives/afni/sub-RG26/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG26	RG	AN	csav	derivatives/afni/sub-RG26/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG26	RG	AN	csunav	derivatives/afni/sub-RG26/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG27	RG	OM	csav	derivatives/afni/sub-RG27/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG27	RG	OM	csm	derivatives/afni/sub-RG27/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG27	RG	AN	csm	derivatives/afni/sub-RG27/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG27	RG	AN	csav	derivatives/afni/sub-RG27/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG27	RG	AN	csunav	derivatives/afni/sub-RG27/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG30	RG	OM	csav	derivatives/afni/sub-RG30/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG30	RG	OM	csm	derivatives/afni/sub-RG30/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG30	RG	AN	csm	derivatives/afni/sub-RG30/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG30	RG	AN	csav	derivatives/afni/sub-RG30/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG30	RG	AN	csunav	derivatives/afni/sub-RG30/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG31	RG	OM	csav	derivatives/afni/sub-RG31/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG31	RG	OM	csm	derivatives/afni/sub-RG31/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG31	RG	AN	csm	derivatives/afni/sub-RG31/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG31	RG	AN	csav	derivatives/afni/sub-RG31/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG31	RG	AN	csunav	derivatives/afni/sub-RG31/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG32	RG	OM	csav	derivatives/afni/sub-RG32/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG32	RG	OM	csm	derivatives/afni/sub-RG32/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG32	RG	AN	csm	derivatives/afni/sub-RG32/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG32	RG	AN	csav	derivatives/afni/sub-RG32/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG32	RG	AN	csunav	derivatives/afni/sub-RG32/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG33	RG	OM	csav	derivatives/afni/sub-RG33/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG33	RG	OM	csm	derivatives/afni/sub-RG33/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG33	RG	AN	csm	derivatives/afni/sub-RG33/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG33	RG	AN	csav	derivatives/afni/sub-RG33/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG33	RG	AN	csunav	derivatives/afni/sub-RG33/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG34	RG	OM	csav	derivatives/afni/sub-RG34/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG34	RG	OM	csm	derivatives/afni/sub-RG34/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG34	RG	AN	csm	derivatives/afni/sub-RG34/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG34	RG	AN	csav	derivatives/afni/sub-RG34/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG34	RG	AN	csunav	derivatives/afni/sub-RG34/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG35	RG	OM	csav	derivatives/afni/sub-RG35/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG35	RG	OM	csm	derivatives/afni/sub-RG35/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG35	RG	AN	csm	derivatives/afni/sub-RG35/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG35	RG	AN	csav	derivatives/afni/sub-RG35/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG35	RG	AN	csunav	derivatives/afni/sub-RG35/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG37	RG	OM	csav	derivatives/afni/sub-RG37/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG37	RG	OM	csm	derivatives/afni/sub-RG37/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG37	RG	AN	csm	derivatives/afni/sub-RG37/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG37	RG	AN	csav	derivatives/afni/sub-RG37/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG37	RG	AN	csunav	derivatives/afni/sub-RG37/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG38	RG	OM	csav	derivatives/afni/sub-RG38/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG38	RG	OM	csm	derivatives/afni/sub-RG38/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG38	RG	AN	csm	derivatives/afni/sub-RG38/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG38	RG	AN	csav	derivatives/afni/sub-RG38/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG38	RG	AN	csunav	derivatives/afni/sub-RG38/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG41	RG	OM	csav	derivatives/afni/sub-RG41/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG41	RG	OM	csm	derivatives/afni/sub-RG41/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG41	RG	AN	csm	derivatives/afni/sub-RG41/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG41	RG	AN	csav	derivatives/afni/sub-RG41/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG41	RG	AN	csunav	derivatives/afni/sub-RG41/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG44	RG	OM	csav	derivatives/afni/sub-RG44/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG44	RG	OM	csm	derivatives/afni/sub-RG44/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG44	RG	AN	csm	derivatives/afni/sub-RG44/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG44	RG	AN	csav	derivatives/afni/sub-RG44/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG44	RG	AN	csunav	derivatives/afni/sub-RG44/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG47	RG	OM	csav	derivatives/afni/sub-RG47/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG47	RG	OM	csm	derivatives/afni/sub-RG47/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG47	RG	AN	csm	derivatives/afni/sub-RG47/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG47	RG	AN	csav	derivatives/afni/sub-RG47/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG47	RG	AN	csunav	derivatives/afni/sub-RG47/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG48	RG	OM	csav	derivatives/afni/sub-RG48/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG48	RG	OM	csm	derivatives/afni/sub-RG48/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG48	RG	AN	csm	derivatives/afni/sub-RG48/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG48	RG	AN	csav	derivatives/afni/sub-RG48/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG48	RG	AN	csunav	derivatives/afni/sub-RG48/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG49	RG	OM	csav	derivatives/afni/sub-RG49/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG49	RG	OM	csm	derivatives/afni/sub-RG49/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG49	RG	AN	csm	derivatives/afni/sub-RG49/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG49	RG	AN	csav	derivatives/afni/sub-RG49/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG49	RG	AN	csunav	derivatives/afni/sub-RG49/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG50	RG	OM	csav	derivatives/afni/sub-RG50/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG50	RG	OM	csm	derivatives/afni/sub-RG50/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG50	RG	AN	csm	derivatives/afni/sub-RG50/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG50	RG	AN	csav	derivatives/afni/sub-RG50/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG50	RG	AN	csunav	derivatives/afni/sub-RG50/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
sub-RG51	RG	OM	csav	derivatives/afni/sub-RG51/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csav#0_Coef]' \
sub-RG51	RG	OM	csm	derivatives/afni/sub-RG51/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[csm#0_Coef]' \
sub-RG51	RG	AN	csm	derivatives/afni/sub-RG51/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csm#0_Coef]' \
sub-RG51	RG	AN	csav	derivatives/afni/sub-RG51/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csav#0_Coef]' \
sub-RG51	RG	AN	csunav	derivatives/afni/sub-RG51/Avoidance/FL_results_WB_pooled/betas_REML+tlrc'[onset_csunav#0_Coef]' \
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
