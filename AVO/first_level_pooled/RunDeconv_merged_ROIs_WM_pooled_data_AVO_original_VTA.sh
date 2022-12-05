# DESCRIPTION: this script allows only to calculate the beta values for each condition of interest for each phase (pooled csm csav for the whole brain AND each ROI). I did this to have the results from the whole brain> For the betas values from each ROI check the output from the other deconv script.

# One time the script is run, you will open AFNI, read the directory and load the results from the underlay option.

#!/bin/bash

set -e

recompute=true
subjects="sub-FG01 sub-FG02 sub-FG03 sub-FG04 sub-FG05 sub-FG06 sub-FG07 sub-FG09 sub-FG10 sub-FG11 sub-FG12 sub-FG14 sub-FG15 sub-FG16 sub-FG17 sub-FG20 sub-FG22 sub-FG23 sub-FG25 sub-FG28R sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-FG33R sub-RG26 sub-RG27 sub-RG30 sub-RG31 sub-RG32 sub-RG33 sub-RG34 sub-RG35 sub-RG37 sub-RG38 sub-RG41 sub-RG44 sub-RG47 sub-RG48 sub-RG49 sub-RG50 sub-RG51 sub-RG55R sub-RG56R sub-RG58R sub-RG59R sub-RG61R sub-RG62R sub-RG63R sub-RG64R"


ROIs="canlab_VTA_bilateral NAcc canlab_NAC_bilateral VmPFC Nacc_shell Nacc_core"

task="Avoidance"

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

if [ ! -e "derivatives/afni/mergedROIs_resam.nii" ]; then  
3dmask_tool \
  -input "derivatives/ROIs/VmPFC_resam.nii" "derivatives/ROIs/NAcc_resam.nii" "derivatives/ROIs/VTA_resam.nii" \
  -prefix "derivatives/afni/mergedROIs_resam.nii" -union
fi
	
	
for subj in $subjects; do
	prefix="derivatives/afni/$subj/${task}"
	result_VW_mergedROIs_prefix="derivatives/afni/$subj/${task}/FL_results_VW_mergedROIs"
	input="$prefix/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz"
	mask="derivatives/afni/mergedROIs_resam.nii"
	all_confounds="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_desc-confounds_timeseries.tsv"
	event_file="$subj/func/${subj}_task-${task}_events.tsv"
	
	mkdir -p $result_VW_mergedROIs_prefix	
	
Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
X <- X[order(X$onset),] # make sure events are ordered by onset time
X <- split(X, X$trial_type) # split by trial_type
prefix <- args[2]

afni_modulated <- function(Z, fname) cat(paste(Z$onset, Z$relief_rating, sep="*"), sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))
afni_unmodulated <- function(Z, fname) cat(Z$onset, sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))
afni_duration <- function(Z, fname) cat(paste(Z$onset, Z$duration, sep=":"), sep=" ", file=file(paste0(prefix, "/", fname, ".1D")))

afni_unmodulated(X[["omission_csm"]], "relief_csm")
afni_unmodulated(X[["omission_csav"]], "omission_csav")

afni_duration(X[["relief_rating"]], "relief_rating")
afni_duration(X[["onset_csm"]], "onset_csm")
afni_duration(X[["onset_csav"]], "onset_csav")
afni_duration(X[["onset_csunav"]], "onset_csunav")

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



	# create highpass regressors (180s)
	1dBport -input $input -band 0 0.005555555555555556 -nozero > "$prefix/highpass.1D"

	# generate X matrix
	3dDeconvolve -input $input \
	    -mask $mask \
	    -nobout \
	    -polort 0 \
	    -local_times \
	    -GOFORIT 3 \
	    -num_stimts 8 \
	    -stim_times 1 "$prefix/omission_csav.1D" 'GAM' \
	    -stim_label 1 csav \
	    -stim_times 2 "$prefix/relief_csm.1D" 'GAM' \
	    -stim_label 2 csm \
	    -stim_times_AM1 3 "$prefix/relief_rating.1D" 'dmBLOCK(1)' \
	    -stim_label 3 relief_rating \
	    -stim_times_AM1 4 "$prefix/onset_csm.1D" 'dmBLOCK(1)' \
	    -stim_label 4 onset_csm \
	    -stim_times_AM1 5 "$prefix/onset_csav.1D" 'dmBLOCK(1)' \
	    -stim_label 5 onset_csav \
	    -stim_times_AM1 6 "$prefix/onset_csunav.1D" 'dmBLOCK(1)' \
	    -stim_label 6 onset_csunav \
	    -stim_times 7 "$prefix/shock.1D" 'GAM' \
	    -stim_label 7 shock \
	    -stim_times 8 "$prefix/press.1D" 'GAM' \
	    -stim_label 8 press \
	    -ortvec $prefix/confounds "confounds" \
	    -ortvec $prefix/highpass.1D highpass \
	    -gltsym 'SYM: csav -csm' \
	    -glt_label 1 "csav vs csm" \
	    -gltsym 'SYM: csav +csm' \
	    -glt_label 2 "csav and csm" \
	    -gltsym 'SYM: csm -csav' \
	    -glt_label 3 "csm vs csav" \
	    -gltsym 'SYM: onset_csav -onset_csm' \
	    -glt_label 4 "onset_csav vs onset_csm" \
	    -gltsym 'SYM: onset_csm -onset_csav' \
	    -glt_label 5 "onset_csm vs onset_csav" \
	    -gltsym 'SYM: onset_csm -onset_csunav' \
	    -glt_label 6 "onset_csm vs onset_csunav" \
	    -gltsym 'SYM: onset_csav -onset_csunav'\
	    -glt_label 7 "onset_csav vs onset_csunav" \
	    -jobs 8 \
	    -allzero_OK \
	    -x1D $result_VW_mergedROIs_prefix/X.xmat.1D \
	    -xjpeg $result_VW_mergedROIs_prefix/X.jpg \
	    -bucket $result_VW_mergedROIs_prefix/betas_ROI \
	    -cbucket $result_VW_mergedROIs_prefix/decon_WB \
	    -tout \
	    -x1D_stop

	[ -e $result_VW_mergedROIs_prefix/betas_REML+tlrc.HEAD -a $recompute == "true" ] && rm "$result_VW_mergedROIs_prefix"/betas_REML*
	if [ ! -e $result_VW_mergedROIs_prefix/betas_REML+tlrc.HEAD ]; then
		3dREMLfit -matrix $result_VW_mergedROIs_prefix/X.xmat.1D \
		     -input $input -mask $mask \
		     -GOFORIT \
		     -Rbuck "$result_VW_mergedROIs_prefix/betas_REML" \
		     -Rvar "$result_VW_mergedROIs_prefix/betas_REMLvar" \
		     -Grid 3 -tout
	fi

	[ -e 3dDeconvolve.err ] && mv 3dDeconvolve.err "$result_VW_mergedROIs_prefix/3dDeconvolve.err"
	[ -e 3dREMLfit.err ] && mv 3dREMLfit.err "$result_VW_mergedROIs_prefix/3dREMLfit.err"
done

