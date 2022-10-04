# DESCRIPTION: this script allows only to calculate the beta values for each condition of interest for each phase (pooled csm csav for the whole brain AND each ROI). I did this to have the results from the whole brain> For the betas values from each ROI check the output from the other deconv script.

# One time the script is run, you will open AFNI, read the directory and load the results from the underlay option.

#!/bin/bash

set -e

recompute=true
subjects="sub-FG01 sub-FG02 sub-FG03 sub-FG04 sub-FG05 sub-FG06 sub-FG07 sub-FG09 sub-FG10 sub-FG11 sub-FG12 sub-FG14 sub-FG15 sub-FG16 sub-FG17 sub-FG20 sub-FG22 sub-FG23 sub-FG25 sub-FG28R sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-FG33R sub-RG26 sub-RG27 sub-RG30 sub-RG31 sub-RG32 sub-RG33 sub-RG34 sub-RG35 sub-RG37 sub-RG38 sub-RG41 sub-RG44 sub-RG47 sub-RG48 sub-RG49 sub-RG50 sub-RG51 sub-RG55R sub-RG56R sub-RG58R sub-RG59R sub-RG61R sub-RG62R sub-RG63R sub-RG64R"

ROIs="VTA NAcc VmPFC VmPFC_box"

task="Avoidance"

wm_inv_mask=derivatives/ROIs/WM_inv_mask.nii.gz

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

Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
X <- X[order(X$onset),] # make sure events are ordered by onset time
phase <- factor((X$trial_number_cs - 1) %/% 4, levels=c(0,1,2), labels=c("early", "middle", "late"))
Y <- split(X, phase)
Y <- lapply(Y, function(x) split(x, x$trial_type))# split by trial_type
X <- split(X, X$trial_type) # split by trial_type
prefix <- args[2]

afni_modulated <- function(Z, fname) cat(paste(Z$onset, Z$relief_rating, sep="*"), sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))
afni_unmodulated <- function(Z, fname) {
  if(nrow(Z) == 0) cat("*", file=file(paste0(prefix, "/", fname, ".1D")))
  else cat(Z$onset, sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))
  
}
afni_duration <- function(Z, fname) cat(paste(Z$onset, Z$duration, sep=":"), sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))


afni_unmodulated(Y[["early"]][["omission_csm"]], "omission_csm_early")
afni_unmodulated(Y[["late"]][["omission_csm"]], "omission_csm_late")
afni_unmodulated(Y[["early"]][["omission_csav"]], "omission_csav_early")
afni_unmodulated(Y[["late"]][["omission_csav"]], "omission_csav_late")

afni_duration(rbind(Y[["early"]][["relief_rating"]],X[["late"]][["relief_rating"]]), "relief_rating")
afni_duration(Y[["early"]][["onset_csm"]], "onset_csm_early")
afni_duration(Y[["late"]][["onset_csm"]], "onset_csm_late")
afni_duration(Y[["early"]][["onset_csav"]], "onset_csav_early")
afni_duration(Y[["late"]][["onset_csav"]], "onset_csav_late")
afni_duration(Y[["early"]][["onset_csunav"]], "onset_csunav_early")
afni_duration(Y[["late"]][["onset_csunav"]], "onset_csunav_late")

afni_unmodulated(rbind(X[["press_csav"]], X[["press_csm"]], X[["press_csunav"]]), "press")
afni_unmodulated(X[["shock"]], "shock")
' $event_file $prefix

# Creates the confounds
Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
n <- names(X)

#CONFOUNDS
S <- startsWith(n, "trans") | startsWith(n, "rot") | n == "csf" | startsWith(n, "motion_outlier")
X <- lapply(X[S], function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
# tried on FG01: these confounds produced activations outside the brain and in the crf,(visible from the f-stat contrast)


write.table(X, file=args[2], row.names=F, col.names=F, sep=" ")
' $all_confounds $prefix/confounds

	# resample the individual gm mask
	if [ ! -e $mask ]; then
		3dresample -master $input \
		  -prefix $mask \
		  -inset $input_mask \
		  -rmode NN
	fi

	# resample WM_inv
	if [ ! -e "$prefix/WM_inv_resam.nii" ]; then
		3dresample -master $input \
		  -prefix $intermediate_tmp \
		  -prefix "$prefix/WM_inv_resam.nii" \
		  -inset $wm_inv_mask \
		  -rmode NN
	fi

# Create ROIs for subject
	if [ ! -e "$prefix/VTA_resam.nii" ]; then
	# resample anat_roi to same resolution as master (your functional images [that you can check it via MANGO, open the func AND ROI, ctrl+I--> image dimension: are they the same??])
		3dresample -master $input \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/VTA_bram.nii.gz" \
		  -rmode NN
		mv $intermediate_tmp "$prefix/VTA_resam.nii" # originql VTA from Esser (VTA plus SN)
	fi

	if [ ! -e "$prefix/NAcc_resam.nii" ]; then
	# resample anat_roi to same resolution as master (your functional images [that you can check it via MANGO, open the func AND ROI, ctrl+I--> image dimension: are they the same??])
		3dresample -master $input \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/NAcc_bil.nii.gz" \
		  -rmode NN

		3dmask_tool -input "$prefix/WM_inv_resam.nii" $intermediate_tmp -prefix "$prefix/NAcc_bil.nii" -inter
		rm $intermediate_tmp
	fi

	if [ ! -e "$prefix/VmPFC_resam.nii" ]; then
	# resample anat_roi to same resolution as master (your functional images [that you can check it via MANGO, open the func AND ROI, ctrl+I--> image dimension: are they the same??])
		3dresample -master $input \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/VmPFC_parcels.nii.gz" \
		  -rmode NN

		3dmask_tool -input "$prefix/WM_inv_resam.nii" $intermediate_tmp -prefix "$prefix/VmPFC_resam.nii" -inter
		rm $intermediate_tmp
		
	fi
		if [ ! -e "$prefix/VmPFC_box_resam.nii" ]; then
	# resample anat_roi to same resolution as master (your functional images [that you can check it via MANGO, open the func AND ROI, ctrl+I--> image dimension: are they the same??])
		3dresample -master $input \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/VmPFC_box.nii" \
		  -rmode NN

		mv $intermediate_tmp "$prefix/VmPFC_box_resam.nii"
		fi

	# create highpass regressors (180s)
	1dBport -input $input -band 0 0.005 -nozero > "$prefix/highpass.1D"

	# generate X matrix
	3dDeconvolve -input $input \
	    -mask $mask \
	    -nobout \
	    -polort 0 \
	    -local_times \
	    -GOFORIT 3 \
	    -num_stimts 13 \
	    -stim_times 1 "$prefix/omission_csav_early.1D" 'BLOCK(4.5)' \
	    -stim_label 1 csav_early \
	    -stim_times 2 "$prefix/omission_csm_early.1D"  'BLOCK(4.5)' \
	    -stim_label 2 csm_early \
	    -stim_times 3 "$prefix/omission_csav_late.1D"  'BLOCK(4.5)' \
	    -stim_label 3 csav_late \
	    -stim_times 4 "$prefix/omission_csm_late.1D"  'BLOCK(4.5)' \
	    -stim_label 4 csm_late \
	    -stim_times_AM1 5 "$prefix/relief_rating.1D" 'dmBLOCK(1)' \
	    -stim_label 5 relief_rating \
	    -stim_times_AM1 6 "$prefix/onset_csm_early.1D" 'dmBLOCK(1)' \
	    -stim_label 6 onset_csm_early \
	    -stim_times_AM1 7 "$prefix/onset_csav_early.1D" 'dmBLOCK(1)' \
	    -stim_label 7 onset_csav_early \
	    -stim_times_AM1 8 "$prefix/onset_csunav_early.1D" 'dmBLOCK(1)' \
	    -stim_label 8 onset_csunav_early \
	    -stim_times_AM1 9 "$prefix/onset_csm_late.1D" 'dmBLOCK(1)' \
	    -stim_label 9 onset_csm_late \
	    -stim_times_AM1 10 "$prefix/onset_csav_late.1D" 'dmBLOCK(1)' \
	    -stim_label 10 onset_csav_late \
	    -stim_times_AM1 11 "$prefix/onset_csunav_late.1D" 'dmBLOCK(1)' \
	    -stim_label 11 onset_csunav_late \
	    -stim_times 12 "$prefix/shock.1D" 'GAM' \
	    -stim_label 12 shock \
	    -stim_times 13 "$prefix/press.1D" 'GAM' \
	    -stim_label 13 press \
	    -ortvec $prefix/confounds "confounds" \
	    -ortvec $prefix/highpass.1D highpass \
	    -gltsym 'SYM: csav_early -csm_early' \
	    -glt_label 1 "csav_early vs csm_early" \
	    -gltsym 'SYM: csav_late -csm_late' \
	    -glt_label 2 "csav_late vs csm_late" \
	    -gltsym 'SYM: csm_early -csav_early' \
	    -glt_label 3 "csm_early vs csav_early" \
	    -gltsym 'SYM: csm_late -csav_late' \
	    -glt_label 4 "csm_late vs csav_late" \
	    -gltsym 'SYM: onset_csav_early -onset_csm_early' \
	    -glt_label 5 "onset_csav_early vs onset_csm_early" \
	    -gltsym 'SYM: onset_csm_early -onset_csav_early' \
	    -glt_label 6 "onset_csm_early vs onset_csav_early" \
	    -gltsym 'SYM: onset_csav_late -onset_csm_late' \
	    -glt_label 7 "onset_csav_late vs onset_csm_late" \
	    -gltsym 'SYM: csav_early -csm_early -csav_late +csm_late' \
	    -glt_label 8 "csav_early_vs_csm_early_minus_csav_late_vs_csm_late" \
	    -gltsym 'SYM: onset_csav_early -onset_csav_late -onset_csm_early +onset_csm_late' \
	    -glt_label 9 "onset_csav_early_vs_late vs onset_csm_early_vs_onset_csm_late" \
	    -gltsym 'SYM: csav_early +csm_early - csav_late -csm_late' \
	    -glt_label 10 "csav_early AND csm_early vs csav_late AND csm_late" \
	    -jobs 8 \
	    -x1D $result_WB_prefix/X.xmat.1D \
	    -xjpeg $result_WB_prefix/X.jpg \
	    -bucket $result_WB_prefix/betas_ROI \
	    -tout \
	    -x1D_stop

	[ -e $result_WB_prefix/betas_REML+tlrc.BRIK -a $recompute == "zzatrue" ] && rm $result_WB_prefix/betas_REML*
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
		roi_mask="${prefix}/${roi}_resam.nii"
		averaged_BOLD_from_ROI="$result_ROI_prefix/averaged_BOLD_from_${roi}.1D"

	        [ -e "$result_ROI_prefix/betas_ROI_${roi}_REML.1D" -a $recompute == "true" ] && rm "$result_ROI_prefix/betas_ROI_${roi}_REML.1D" "$result_ROI_prefix/betas_ROI_${roi}_REMLvar.1D"
		if [ ! -e "$result_ROI_prefix/betas_ROI_${roi}_REML.1D" ]; then
			# average voxel signal to get the mean betas using ROIs before using the events
			3dmaskave -quiet -mask $roi_mask $input > $averaged_BOLD_from_ROI

			3dREMLfit -matrix $result_WB_prefix/X.xmat.1D \
			     -input ${averaged_BOLD_from_ROI}'[0]'\' \
			     -Rbuck "$result_ROI_prefix/betas_ROI_${roi}_REML" \
			     -Rvar "$result_ROI_prefix/betas_ROI_${roi}_REMLvar" \
			     -GOFORIT \
			     -Grid 3 -tout
		fi
	done
done


