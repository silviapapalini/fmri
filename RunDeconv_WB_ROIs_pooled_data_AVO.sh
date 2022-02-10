# DESCRIPTION: this script allows only to calculate the beta values for each condition of interest for each phase (pooled csm csav for the whole brain AND each ROI). I did this to have the results from the whole brain> For the betas values from each ROI check the output from the other deconv script.

# One time the script is run, you will open AFNI, read the directory and load the results from the underlay option.

#!/bin/bash

set -e

#subjects="sub-FG01"
subjects="sub-FG01 sub-FG02 sub-FG03 sub-FG04 sub-FG05 sub-FG06 sub-FG07 sub-FG09 sub-FG10 sub-FG11 sub-FG12 sub-FG14 sub-FG15 sub-FG16 sub-FG17 sub-FG20 sub-FG22 sub-FG23 sub-FG25 sub-FG28R sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-FG33R sub-RG26 sub-RG27 sub-RG30 sub-RG31 sub-RG32 sub-RG33 sub-RG34 sub-RG35 sub-RG37 sub-RG38 sub-RG41 sub-RG44 sub-RG47 sub-RG48 sub-RG49 sub-RG50 sub-RG51 sub-RG55R sub-RG56R sub-RG58R sub-RG59R sub-RG61R sub-RG62R sub-RG63R sub-RG64R"

#subjects="sub-RG41 sub-RG44 sub-RG47 sub-RG48 sub-RG49 sub-RG50 sub-RG51 sub-RG55R sub-RG56R sub-RG58R sub-RG59R sub-RG61R sub-RG62R sub-RG63R sub-RG64R"


ROIs="VTA NAcc VmPFC"

task="Avoidance"
# resample anat_roi to same resolution as master (your functional images [that you can check it via MANGO, open the func AND ROI, ctrl+I--> image dimension: are they the same??])

input=derivatives/afni/sub-FG01/Extinction/sub-FG01_task-Extinction_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz 
if [ ! -e "derivatives/ROIs/VTA_resam.nii" ]; then
	3dresample -master $input \
	  -prefix "derivatives/ROIs/VTA_resam.nii" \
	  -inset "derivatives/ROIs/SNVTA_MTmean.hdr" \
	  -rmode NN
fi

if [ ! -e "derivatives/ROIs/NAcc_resam.nii" ]; then
	3dresample -master $input \
	  -prefix "derivatives/ROIs/NAcc_resam.nii" \
	  -inset "derivatives/ROIs/NAcc_dilated.nii" \
	  -rmode NN
fi

if [ ! -e "derivatives/ROIs/VmPFC_resam.nii" ]; then  
	3dresample -master $input \
	  -prefix "derivatives/ROIs/VmPFC_resam.nii" \
	  -inset "derivatives/ROIs/VmPFC.nii" \
	  -rmode NN
fi

for subj in $subjects; do
	prefix="derivatives/afni/$subj/${task}"
	result_WB_prefix="derivatives/afni/$subj/${task}/FL_results_WB_pooled"
	result_ROI_prefix="derivatives/afni/$subj/${task}/FL_results_ROI_pooled"
	input="$prefix/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz"
	mask="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"
	all_confounds="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_desc-confounds_timeseries.tsv"
	event_file="$subj/func/${subj}_task-${task}_events.tsv"
	
	mkdir -p $result_WB_prefix $result_ROI_prefix 
	
Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
X <- X[order(X$onset),] # make sure events are ordered by onset time
X <- split(X, X$trial_type) # split by trial_type
prefix <- args[2]
cat(X[["omission_csm_onset_relief"]]$onset, sep=" ", file=file(paste0(prefix, "/relief_csm.1D")))
cat(X[["omission_csav_onset_relief"]]$onset, sep=" ", file=file(paste0(prefix, "/relief_csav.1D")))
cat(X[["onset_rating_relief"]]$onset, sep=" ", file=file(paste0(prefix, "/relief_rating.1D")))
' $event_file $prefix
	
	Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
n <- names(X)


#I CAN CHOOSE between:

#CONFOUNDS 1
S <- startsWith(n, "trans") | startsWith(n, "rot") | startsWith(n, "white_matter") | startsWith(n, "csf") | startsWith(n, "motion_outlier")
#I did not include Global signal, and so there are 33 regressors for movment. I tried also | startsWith (n, "cosine00"), but results are also worse 

#CONFOUNDS 2
#S <- startsWith(n, "trans") | startsWith(n, "rot") | n == "white_matter" | n == "csf" | startsWith(n, "motion_outlier")
X <- lapply(X[S], function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
# tried on FG01: these confounds produced activations outside the brain and in the crf,(visible from the f-stat contrast)


write.table(X, file=args[2], row.names=F, col.names=F, sep=" ")
' $all_confounds $prefix/confounds
	
	if [ ! -e $result_WB_prefix/statsWB+tlrc.BRIK ]; then
		3dDeconvolve -input $input \
			-mask $mask \
			-polort A \
			-goforit 2 \
			-local_times \
			-num_stimts 3 \
			-stim_times 1 "$prefix/relief_csav.1D" 'GAM' \
			-stim_label 1 csav \
			-stim_times 2 "$prefix/relief_csm.1D" 'GAM' \
			-stim_label 2 csm \
			-stim_times 3 "$prefix/relief_rating.1D" 'GAM' \
			-stim_label 3 relief_rating \
			-ortvec $prefix/confounds "confounds" \
			-jobs 8 \
			-gltsym 'SYM: csav -csm' \
			-glt_label 1 csav-csm  \
			-gltsym 'SYM: csm -csav' \
			-glt_label 2 csm-csav \
			-gltsym 'SYM: csav' \
			-glt_label 3 csav  \
			-gltsym 'SYM: csm' \
			-glt_label 4 csm \
			-fout -tout \
			-x1D $result_WB_prefix/XWB.xmat.1D \
			-xjpeg $result_WB_prefix/XWB.jpg \
			-fitts $result_WB_prefix/fittsWB \
			-errts $result_WB_prefix/errtsWB \
			-bucket $result_WB_prefix/statsWB
		  [ -e 3dDeconvolve.err ] && mv 3dDeconvolve.err "$result_WB_prefix/3dDeconvolveWB.err"
	else
		echo "$result_WB_prefix/statsWB already exists"
	fi
		  
	for roi in $ROIs; do
		roi_mask="derivatives/ROIs/${roi}_resam.nii"
		3dmaskave -quiet -mask $roi_mask $result_WB_prefix/statsWB+tlrc[10,13,16,19] > "$result_ROI_prefix/averaged_betas_from_${roi}.1D"
	done
done


