# READ ME
# This script allows only to calculate the voxelwise analysis for each condition and trial from the VmPFC (csm csav right and left unified). 
# Since the event files are splitted there is not real info in the output matrix about the order of presentation of the events. This means that before to run the LMM you will need to provide the real order of presentation of each CS as it was in the task/learning phase.
# One time the script is run, you will open AFNI, read the directory and load the betas from the underlay option.

#!/bin/bash

set -e

recompute=true
subjects="sub-FG01 sub-FG02 sub-FG03 sub-FG04 sub-FG05 sub-FG06 sub-FG07 sub-FG09 sub-FG10 sub-FG11 sub-FG12 sub-FG14 sub-FG16 sub-FG17 sub-FG20 sub-FG22 sub-FG23 sub-FG25 sub-FG28R sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-FG33R sub-RG26 sub-RG27 sub-RG30 sub-RG31 sub-RG34 sub-RG37 sub-RG38 sub-RG41 sub-RG44 sub-RG47 sub-RG48 sub-RG49 sub-RG50 sub-RG51 sub-RG55R sub-RG56R sub-RG58R sub-RG59R sub-RG61R sub-RG62R sub-RG63R sub-RG64R"

# select/deselect the learning phase of interest. PAV, AVO, EXT
task="Extinction"

# resample anat_roi to same resolution as master (your functional images [that you can check it via MANGO, open the func AND ROI, ctrl+I--> image dimension: are they the same??])
input=derivatives/afni/sub-FG01/Extinction/sub-FG01_task-Extinction_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz 

if [ ! -e "derivatives/ROIs/VmPFC_resam.nii" ]; then  
3dresample -master $input \
  -prefix "derivatives/ROIs/VmPFC_resam.nii" \
  -inset "derivatives/ROIs/VmPFC_parcels.nii.gz" \
  -rmode NN
fi


for subj in $subjects; do
	prefix="derivatives/afni/$subj/${task}"
	result_VW_VmPFC_prefix="derivatives/afni/$subj/${task}/FL_results_VW_VmPFC_parametric_modulation"
	input="$prefix/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz"
	mask="derivatives/ROIs/VmPFC_resam.nii"
	all_confounds="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_desc-confounds_timeseries.tsv"
	event_file="$subj/func/${subj}_task-${task}_events.tsv"
	
	mkdir -p $result_VW_VmPFC_prefix

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
            -gltsym 'SYM: csm_modulated[0] -csav_modulated[0]' \
            -glt_label 4 "csmnonM - csavnonM" \
	    -gltsym 'SYM: csm_modulated[1] -csav_modulated[1]'  \
            -glt_label 5 "csm_modulated - csav_modulated" \
            -gltsym 'SYM: csm_modulated[1] +csav_modulated[1]'  \
            -glt_label 6 "csm_modulated AND csav_modulated" \
	    -jobs 8 \
	    -allzero_OK \
	    -x1D $result_VW_VmPFC_prefix/X.xmat.1D \
	    -xjpeg $result_VW_VmPFC_prefix/X.jpg \
	    -bucket $result_VW_VmPFC_prefix/betas_ROI \
	    -cbucket $result_VW_VmPFC_prefix/decon_WB \
	    -tout \
	    -x1D_stop

	[ -e $result_VW_VmPFC_prefix/betas_REML+tlrc.HEAD -a $recompute == "true" ] && rm "$result_VW_VmPFC_prefix"/betas_REML*
	if [ ! -e $result_VW_VmPFC_prefix/betas_REML+tlrc.HEAD ]; then
		3dREMLfit -matrix $result_VW_VmPFC_prefix/X.xmat.1D \
		     -input $input -mask $mask \
		     -GOFORIT \
		     -Rbuck "$result_VW_VmPFC_prefix/betas_REML" \
		     -Rvar "$result_VW_VmPFC_prefix/betas_REMLvar" \
		     -Grid 3 -tout
	fi

	[ -e 3dDeconvolve.err ] && mv 3dDeconvolve.err "$result_VW_VmPFC_prefix/3dDeconvolve.err"
	[ -e 3dREMLfit.err ] && mv 3dREMLfit.err "$result_VW_VmPFC_prefix/3dREMLfit.err"
done

