# READ ME
# This script allows only to calculate the beta values for each condition and trial from each ROI (csm csav for the VTA, Nacc, vmPFC, right and left unified). 
# Since the event files are splitted there is not real info in the output matrix about the order of presentation of the events. This means that before to run the LMM you will need to provide the real order of presentation of each CS as it was in the task/learning phase.
# One time the script is run, you will open AFNI, read the directory and load the betas from the underlay option.

#!/bin/bash

set -e

recompute=true
subjects="sub-FG01 sub-FG02 sub-FG03 sub-FG04 sub-FG05 sub-FG06 sub-FG07 sub-FG09 sub-FG10 sub-FG11 sub-FG12 sub-FG14 sub-FG16 sub-FG17 sub-FG20 sub-FG22 sub-FG23 sub-FG25 sub-FG28R sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-FG33R sub-RG26 sub-RG27 sub-RG30 sub-RG31 sub-RG34 sub-RG37 sub-RG38 sub-RG41 sub-RG44 sub-RG47 sub-RG48 sub-RG49 sub-RG50 sub-RG51 sub-RG55R sub-RG56R sub-RG58R sub-RG59R sub-RG61R sub-RG62R sub-RG63R sub-RG64R"

# select/deselect the learning phase of interest. PAV, AVO, EXT
task="Extinction"
ROIs="VTA NAcc VmPFC"

wm_inv_mask=derivatives/ROIs/WM_inv_mask.nii.gz

intermediate_tmp=$(mktemp -u --suffix=.nii)
for subj in $subjects; do
	prefix="derivatives/afni/$subj/${task}"
	result_WB_prefix="derivatives/afni/$subj/${task}/FL_results_WB_parametric_modulation"
	result_ROI_prefix="derivatives/afni/$subj/${task}/FL_results_ROI_WM_parametric_modulation"
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
X <- split(X, X$trial_type) # split by trial_type
prefix <- args[2]

afni_modulated <- function(Z, fname) cat(paste(Z$onset, Z$relief_rating, sep="*"), sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))
afni_unmodulated <- function(Z, fname) cat(Z$onset, sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))
afni_duration <- function(Z, fname) cat(paste(Z$onset, Z$duration, sep=":"), sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))

afni_modulated(X[["relief_csm"]], "relief_csm_modulated")
afni_modulated(X[["relief_csav"]], "relief_csav_modulated")

afni_duration(X[["relief_rating"]], "relief_rating")
afni_duration(X[["onset_csm"]], "onset_csm")
afni_duration(X[["onset_csav"]], "onset_csav")
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

	# resample WM_inv
	if [ ! -e "$prefix/WM_inv_resam.nii" ]; then
		3dresample -master $input \
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
		  -inset "derivatives/ROIs/NAcc_HarvardOxford.nii.gz" \
		  -rmode NN

		3dmask_tool -input "$prefix/WM_inv_resam.nii" $intermediate_tmp -prefix "$prefix/NAcc_resam.nii" -inter
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

	# create highpass regressors (180s)
	1dBport -input $input -band 0 0.005 -nozero > "$prefix/highpass.1D"

	# generate X matrix
	3dDeconvolve \
	    -input $input \
	    -mask $mask \
	    -nobout \
	    -polort 0 \
	    -local_times \
	    -num_stimts 5 \
	    -stim_times_AM2 1 "$prefix/relief_csav_modulated.1D" 'BLOCK(4.5)' \
	    -stim_label 1 csav_modulated \
	    -stim_times_AM2 2 "$prefix/relief_csm_modulated.1D" 'BLOCK(4.5)' \
	    -stim_label 2 csm_modulated \
	    -stim_times_AM1 3 "$prefix/relief_rating.1D" 'dmBLOCK(1)' \
	    -stim_label 3 relief_rating \
	    -stim_times_AM1 4 "$prefix/onset_csm.1D" 'dmBLOCK(1)' \
	    -stim_label 4 onset_csm \
	    -stim_times_AM1 5 "$prefix/onset_csav.1D" 'dmBLOCK(1)' \
	    -stim_label 5 onset_csav \
	    -ortvec $prefix/confounds "confounds" \
	    -ortvec $prefix/highpass.1D highpass \
	    -gltsym 'SYM: csav_modulated[0] -csav_modulated[1]'  \
	    -glt_label 1 "csavnonM - csavM" \
	    -gltsym 'SYM: csm_modulated[0] -csm_modulated[1]'  \
	    -glt_label 2 "csmnonM - csmM" \
            -gltsym 'SYM: csav_modulated[0] +csm_modulated[0] -csav_modulated[1] -csm_modulated[1]' \
            -glt_label 3 "csavnonM AND csmnonM - csavM AND csmM"\
            -gltsym 'SYM: csav_modulated[0] -csm_modulated[0]' \
            -glt_label 4 "csavnonM - csmnonM" \
	    -gltsym 'SYM: csm_modulated[1] -csav_modulated[1]'  \
            -glt_label 5 "csm_modulated - csav_modulated" \
            -gltsym 'SYM: csav_modulated[1] -csm_modulated[1]'  \
            -glt_label 6 "csav_modulated - csm_modulated" \
            -gltsym 'SYM: csav_modulated[1] +csm_modulated[1]'  \
            -glt_label 7 "csav_modulated AND csm_modulated" \
	    -jobs 8 \
	    -allzero_OK \
	    -x1D $result_WB_prefix/X.xmat.1D \
	    -xjpeg $result_WB_prefix/X.jpg \
	    -bucket $result_WB_prefix/betas_ROI \
	    -cbucket $result_WB_prefix/decon_WB \
	    -tout \
	    -x1D_stop

	[ -e $result_WB_prefix/betas_REML+tlrc.BRIK -a $recompute == "true" ] && rm "$result_WB_prefix"/betas_REML*
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
			     -GOFORIT \
			     -Rbuck "$result_ROI_prefix/betas_ROI_${roi}_REML" \
			     -Rvar "$result_ROI_prefix/betas_ROI_${roi}_REMLvar" \
			     -Grid 3 -tout
		fi
	done
done

