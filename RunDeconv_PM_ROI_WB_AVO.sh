# READ ME
# This script allows only to calculate the beta values for each condition and trial from each ROI (csm csav for the VTA, Nacc, vmPFC, right and left unified). 
# Since the event files are splitted there is not real info in the output matrix about the order of presentation of the events. This means that before to run the LMM you will need to provide the real order of presentation of each CS as it was in the task/learning phase.
# One time the script is run, you will open AFNI, read the directory and load the betas from the underlay option.

#!/bin/bash

set -e

subjects="sub-FG01 sub-FG02 sub-FG03 sub-FG04 sub-FG05 sub-FG06 sub-FG07 sub-FG09 sub-FG10 sub-FG11 sub-FG12 sub-FG14 sub-FG16 sub-FG17 sub-FG20 sub-FG22 sub-FG23 sub-FG25 sub-FG28R sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-FG33R sub-RG26 sub-RG27 sub-RG30 sub-RG31 sub-RG34 sub-RG37 sub-RG38 sub-RG41 sub-RG44 sub-RG47 sub-RG48 sub-RG49 sub-RG50 sub-RG51 sub-RG55R sub-RG56R sub-RG58R sub-RG59R sub-RG61R sub-RG62R sub-RG63R sub-RG64R"
subjects="sub-FG01"

# select/deselect the learning phase of interest. PAV, AVO, EXT
task="Avoidance"
ROIs="VTA NAcc VmPFC"

# resample anat_roi to same resolution as master (your functional images [that you can check it via MANGO, open the func AND ROI, ctrl+I--> image dimension: are they the same??])
input=derivatives/afni/sub-FG01/Avoidance/sub-FG01_task-Avoidance_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz 
if [ ! -e "derivatives/ROIs/VTA_resam.nii" ]; then
	3dresample -master $input \
	  -prefix "derivatives/ROIs/VTA_resam.nii" \
	  -inset "derivatives/ROIs/VTA_bram.nii.gz" \
	  -rmode NN
fi

if [ ! -e "derivatives/ROIs/NAcc_resam.nii" ]; then
3dresample -master $input \
  -prefix "derivatives/ROIs/NAcc_resam.nii" \
  -inset "derivatives/ROIs/NAcc_HarvardOxford.nii.gz" \
  -rmode NN
fi

if [ ! -e "derivatives/ROIs/VmPFC_resam.nii" ]; then  
3dresample -master $input \
  -prefix "derivatives/ROIs/VmPFC_resam.nii" \
  -inset "derivatives/ROIs/VmPFC_parcels.nii.gz" \
  -rmode NN
fi
  
for subj in $subjects; do
	prefix="derivatives/afni/$subj/${task}"
	result_WB_prefix="derivatives/afni/$subj/${task}/FL_results_WB_parametric_modulation"
	result_ROI_prefix="derivatives/afni/$subj/${task}/FL_results_ROI_parametric_modulation"
	input="$prefix/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz"
	mask="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"
	all_confounds="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_desc-confounds_timeseries.tsv"
	event_file="$subj/func/${subj}_task-${task}_events.tsv"
	
	mkdir -p $result_WB_prefix $result_ROI_prefix 

# Creates the events files (not pooled as it is for fmriprep)	
Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
X <- X[order(X$onset),] # make sure events are ordered by onset time
X <- split(X, X$trial_type) # split by trial_type
prefix <- args[2]

afni_modulated <- function(Z, fname) cat(paste(Z$onset, Z$relief_rating, sep="*"), sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))
afni_unmodulated <- function(Z, fname) cat(Z$onset, sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))

afni_modulated(rbind(X[["relief_csm"]], X[["omission_csm"]]), "relief_csm_modulated")
afni_modulated(X[["omission_csav"]], "omission_csav_modulated")

cat(paste(X[["relief_rating"]]$onset, X[["relief_rating"]]$duration, sep=":"), sep=" ", file=file(paste0(prefix, "/relief_rating.1D")))
cat(paste(X[["onset_csm"]]$onset, X[["onset_csm"]]$duration, sep=":"), sep=" ", file=file(paste0(prefix, "/onset_csm.1D")))
cat(paste(X[["onset_csav"]]$onset, X[["onset_csav"]]$duration, sep=":"), sep=" ", file=file(paste0(prefix, "/onset_csav.1D")))

afni_unmodulated(rbind(X[["press_csav"]], X[["press_csm"]], X[["press_csunav"]]), "press")
afni_unmodulated(X[["shock"]], "shock")
' $event_file $prefix

# Creates the confounds	
Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
n <- names(X)

#CONFOUNDS 
S <- startsWith(n, "trans") | startsWith(n, "rot") | startsWith(n, "white_matter") | startsWith(n, "csf") | startsWith(n, "motion_outlier") |  startsWith(n, "cosine")
X <- lapply(X[S], function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })

write.table(X, file=args[2], row.names=F, col.names=F, sep=" ")
' $all_confounds $prefix/confounds

	# generate X matrix
	3dDeconvolve -input $input \
	    -mask $mask \
	    -nobout \
	    -polort A \
	    -local_times \
	    -num_stimts 7 \
	    -stim_times_AM2 1 "$prefix/omission_csav_modulated.1D" 'GAM' \
	    -stim_label 1 csav_modulated \
	    -stim_times_AM2 2 "$prefix/relief_csm_modulated.1D" 'GAM' \
	    -stim_label 2 csm_modulated \
	    -stim_times_AM1 3 "$prefix/relief_rating.1D" 'dmBLOCK(1)' \
	    -stim_label 3 relief_rating \
	    -stim_times_AM1 4 "$prefix/onset_csm.1D" 'dmBLOCK(1)' \
	    -stim_label 4 onset_csm \
	    -stim_times_AM1 5 "$prefix/onset_csav.1D" 'dmBLOCK(1)' \
	    -stim_label 5 onset_csav \
	    -stim_times 6 "$prefix/shock.1D" 'GAM' \
	    -stim_label 6 shock \
	    -stim_times 7 "$prefix/press.1D" 'GAM' \
	    -stim_label 7 press \
	    -ortvec $prefix/confounds "confounds" \
	    -jobs 8 \
	    -allzero_OK \
	    -x1D $result_WB_prefix/X.xmat.1D \
	    -xjpeg $result_WB_prefix/X.jpg \
	    -bucket $result_WB_prefix/betas_ROI \
	    -x1D_stop

	[ -e $result_WB_prefix/betas_REML+tlrc.BRIK ] && rm $result_WB_prefix/betas_REML+tlrc.*
	if [ ! -e $result_WB_prefix/betas_REML+tlrc.BRIK ]; then
		3dREMLfit -matrix $result_WB_prefix/X.xmat.1D \
		     -input $input -mask $mask \
		     -GOFORIT \
		     -Rbuck "$result_WB_prefix/betas_REML" \
		     -Grid 3 -tout \
		     -gltsym 'SYM: csav_modulated[0] -csav_modulated[1]' "csavnonM - csavM" \
		     -gltsym 'SYM: csav_modulated[0] -csm_modulated[0]' "csavnonM - csmnonM"
	fi

	[ -e 3dDeconvolve.err ] && mv 3dDeconvolve.err "$result_WB_prefix/3dDeconvolve.err"
	[ -e 3dREMLfit.err ] && mv 3dREMLfit.err "$result_WB_prefix/3dREMLfit.err"

	for roi in $ROIs; do
		roi_mask="derivatives/ROIs/${roi}_resam.nii"
		averaged_BOLD_from_ROI="$result_ROI_prefix/averaged_BOLD_from_${roi}.1D"
	
		# average voxel signal to get the mean betas using ROIs before using the events
		3dmaskave -quiet -mask $roi_mask $input > $averaged_BOLD_from_ROI
	  
	        [ -e "$result_ROI_prefix/betas_ROI_${roi}_REML.1D" ] && rm "$result_ROI_prefix/betas_ROI_${roi}_REML.1D"
		3dREMLfit -matrix $result_WB_prefix/X.xmat.1D \
		     -input ${averaged_BOLD_from_ROI}'[0]'\' \
		     -GOFORIT \
		     -Rbuck "$result_ROI_prefix/betas_ROI_${roi}_REML" \
		     -Grid 3 -tout
	done
done




