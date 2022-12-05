# DESCRIPTION: this script allows only to calculate the beta values for each condition of interest for each phase (pooled csm csav for the whole brain AND each ROI). I did this to have the results from the whole brain> For the betas values from each ROI check the output from the other deconv script.

# One time the script is run, you will open AFNI, read the directory and load the results from the underlay option.

#!/bin/bash

set -e

recompute=true
subjects="sub-FG01 sub-FG02 sub-FG03 sub-FG04 sub-FG05 sub-FG06 sub-FG07 sub-FG09 sub-FG10 sub-FG11 sub-FG12 sub-FG14 sub-FG16 sub-FG17 sub-FG20 sub-FG22 sub-FG23 sub-FG25 sub-FG28R sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-FG33R sub-RG26 sub-RG27 sub-RG30 sub-RG31 sub-RG34 sub-RG37 sub-RG38 sub-RG41 sub-RG44 sub-RG47 sub-RG48 sub-RG49 sub-RG50 sub-RG51 sub-RG55R sub-RG56R sub-RG58R sub-RG59R sub-RG61R sub-RG62R sub-RG63R sub-RG64R"


ROIs="VTA canlab_VTA_bilateral NAcc canlab_NAC_bilateral VmPFC VmPFC_box VmPFC_csp_vs_csm"
task="Extinction"

intermediate_tmp=$(mktemp -u --suffix=.nii)
for subj in $subjects; do
	prefix="derivatives/afni/$subj/${task}"
	result_WB_prefix="derivatives/afni/$subj/${task}/FL_results_WB_pooled_early_late"
	result_ROI_prefix="derivatives/afni/$subj/${task}/FL_results_ROI_WM_pooled_early_late"
	input="$prefix/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz"
	input_mask="derivatives/afni/${subj}/gm_vta_mask.nii.gz"
	mask="$prefix/gm_vta_mask_resam.nii.gz"
	all_confounds="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_desc-confounds_timeseries.tsv"
	event_file="$subj/func/${subj}_task-${task}_events.tsv"
	
	mkdir -p $result_WB_prefix $result_ROI_prefix 

# Creates the events files (not pooled as it is for fmriprep)	
Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
X <- X[order(X$onset),] # make sure events are ordered by onset time
X <- split(X, ifelse(X$trial_number_cs <= 4, "early", "late"))
X <- lapply(X, function(x) split(x, x$trial_type))# split by trial_type
prefix <- args[2]

afni_modulated <- function(Z, fname) cat(paste(Z$onset, Z$relief_rating, sep="*"), sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))
afni_unmodulated <- function(Z, fname) cat(Z$onset, sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))
afni_duration <- function(Z, fname) cat(paste(Z$onset, Z$duration, sep=":"), sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))

afni_unmodulated(X[["early"]][["relief_csm"]], "relief_csm_early")
afni_unmodulated(X[["late"]][["relief_csm"]], "relief_csm_late")
afni_unmodulated(X[["early"]][["relief_csav"]], "relief_csav_early")
afni_unmodulated(X[["late"]][["relief_csav"]], "relief_csav_late")

afni_duration(rbind(X[["early"]][["relief_rating"]],X[["late"]][["relief_rating"]]), "relief_rating")
afni_duration(X[["early"]][["onset_csm"]], "onset_csm_early")
afni_duration(X[["early"]][["onset_csav"]], "onset_csav_early")
afni_duration(X[["late"]][["onset_csm"]], "onset_csm_late")
afni_duration(X[["late"]][["onset_csav"]], "onset_csav_late")
' $event_file $prefix

# Creates the confounds	
Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
n <- names(X)

#CONFOUNDS 
S <- startsWith(n, "trans") | startsWith(n, "rot") | n == "csf" | startsWith(n, "motion_outlier")
X <- lapply(X[S], function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })

write.table(X, file=args[2], row.names=F, col.names=F, sep=" ")
' $all_confounds $prefix/confounds

	# resample the individual gm mask
	if [ ! -e $mask ]; then  
		3dresample -master $input \
		  -prefix $mask \
		  -inset $input_mask \
		  -rmode NN
	fi

	# create highpass regressors (180s)
	1dBport -input $input -band 0 0.005555555555555556 -nozero > "$prefix/highpass.1D"

	# generate X matrix
	3dDeconvolve \
	    -input $input \
	    -mask $mask \
	    -nobout \
	    -polort 0 \
	    -local_times \
	    -num_stimts 9 \
	    -num_glt 13 \
	    -stim_times 1 "$prefix/relief_csav_early.1D" 'BLOCK(4.5)' \
	    -stim_label 1 csav_early \
	    -stim_times 2 "$prefix/relief_csm_early.1D" 'BLOCK(4.5)' \
	    -stim_label 2 csm_early \
	    -stim_times 3 "$prefix/relief_csav_late.1D" 'BLOCK(4.5)' \
	    -stim_label 3 csav_late \
	    -stim_times 4 "$prefix/relief_csm_late.1D" 'BLOCK(4.5)' \
	    -stim_label 4 csm_late \
	    -stim_times_AM1 5 "$prefix/relief_rating.1D" 'dmBLOCK(1)' \
	    -stim_label 5 relief_rating \
	    -stim_times_AM1 6 "$prefix/onset_csm_early.1D" 'dmBLOCK(1)' \
	    -stim_label 6 onset_csm_early \
	    -stim_times_AM1 7 "$prefix/onset_csav_early.1D" 'dmBLOCK(1)' \
	    -stim_label 7 onset_csav_early \
	    -stim_times_AM1 8 "$prefix/onset_csm_late.1D" 'dmBLOCK(1)' \
	    -stim_label 8 onset_csm_late \
	    -stim_times_AM1 9 "$prefix/onset_csav_late.1D" 'dmBLOCK(1)' \
	    -stim_label 9 onset_csav_late \
	    -ortvec $prefix/confounds "confounds" \
	    -ortvec $prefix/highpass.1D highpass \
	    -gltsym 'SYM: csav_early -csm_early ' \
	    -glt_label 1 "csav_early  vs csm_early " \
	    -gltsym 'SYM: csm_early  -csav_early ' \
	    -glt_label 2 "csm_early  vs csav_early " \
	    -gltsym 'SYM: csav_late -csm_late ' \
	    -glt_label 3 "csav_late  vs csm_late " \
	    -gltsym 'SYM: csm_late  -csav_late ' \
	    -glt_label 4 "csm_late  vs csav_late " \
	    -gltsym 'SYM: onset_csav_early -onset_csm_early' \
	    -glt_label 5 "onset_csav_early vs onset_csm_early" \
	    -gltsym 'SYM: onset_csm_early -onset_csav_early' \
	    -glt_label 6 "onset_csm_early vs onset_csav_early" \
            -gltsym 'SYM: onset_csav_late -onset_csm_late' \
	    -glt_label 7 "onset_csav_late vs onset_csm_late" \
	    -gltsym 'SYM: onset_csm_late -onset_csav_late' \
	    -glt_label 8 "onset_csm_late vs onset_csav_late" \
	    -gltsym 'SYM: onset_csav_early -onset_csav_late -onset_csm_early +onset_csm_late' \
	    -glt_label 9 "onset_csav_early_vs_late vs onset_csm_early_vs_onset_csm_late" \
	    -gltsym 'SYM: csav_early +csm_early -csav_late -csm_late' \
	    -glt_label 10 "csav_early_AND_csm_early vs csav_late_AND_csm_late" \
	    -gltsym 'SYM: onset_csav_early +onset_csm_early -onset_csav_late -onset_csm_late' \
	    -glt_label 11 "onset_csav_early_AND_onset_csm_early vs onset_csav_late_AND_onset_csm_late" \
	    -gltsym 'SYM: csm_early -csav_early -csm_late +csav_late' \
	    -glt_label 12 "csm_early_vs_csav_early_minus_csm_late_vs_csav_late" \
	    -gltsym 'SYM: csav_early -csm_early -csav_late +csm_late' \
	    -glt_label 13 "csav_early_vs_csm_early_minus_csav_late_vs_csm_late" \
	    -jobs 8 \
	    -tout \
	    -x1D $result_WB_prefix/X.xmat.1D \
	    -xjpeg $result_WB_prefix/X.jpg \
	    -fitts $result_WB_prefix/fittsWB \
	    -errts $result_WB_prefix/errtsWB \
	    -bucket $result_WB_prefix/statsWB \
	    -x1D_stop

	[ -e $result_WB_prefix/betas_REML+tlrc.BRIK -a $recompute == "zztrue" ] && rm "$result_WB_prefix"/betas_REML*
	if [ ! -e $result_WB_prefix/betas_REML+tlrc.BRIK ]; then
		3dREMLfit -matrix $result_WB_prefix/X.xmat.1D \
		     -input $input -mask $mask \
		     -GOFORIT \
		     -Rbuck "$result_WB_prefix/betas_REML" \
		     -Rvar "$result_WB_prefix/betas_REMLvar" \
		     -Grid 3 -tout
	fi

	[ -e 3dDeconvolve.err ] && mv 3dDeconvolve.err "$result_WB_prefix/3dDeconvolve.err"
	[ -e 3dREMLfit.err ] && mv 3dREMLfit.err "$result_WB_prefix/3dREMLfit.err"

	for roi in $ROIs; do
		#averaged_BOLD_from_ROI="$result_ROI_prefix/averaged_BOLD_from_${roi}.1D"
	        scaled_BOLD_from_ROI="$prefix/${subj}_task-${task}_scaled_BOLD_from_${roi}.1D"
	        [ -e "$result_ROI_prefix/betas_ROI_${roi}_REML.1D" -a $recompute == "true" ] && rm "$result_ROI_prefix/betas_ROI_${roi}_REML.1D" "$result_ROI_prefix/betas_ROI_${roi}_REMLvar.1D"
		if [ ! -e "$result_ROI_prefix/betas_ROI_${roi}_REML.1D" ]; then
			# average voxel signal to get the mean betas using ROIs before using the events

			3dREMLfit -matrix $result_WB_prefix/X.xmat.1D \
			     -input ${scaled_BOLD_from_ROI} \
			     -GOFORIT \
			     -Rbuck "$result_ROI_prefix/betas_ROI_${roi}_REML" \
			     -Rvar "$result_ROI_prefix/betas_ROI_${roi}_REMLvar" \
			     -Grid 3 -tout
		fi
	done

done


