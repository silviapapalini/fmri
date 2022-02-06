# READ ME
# This script allows only to calculate the beta values for each condition and trial from each ROI (csm csav for the VTA, Nacc, vmPFC, right and left unified). 
# Since the event files are splitted there is not real info in the output matrix about the order of presentation of the events. This means that before to run the LMM you will need to provide the real order of presentation of each CS as it was in the task/learning phase.
# One time the script is run, you will open AFNI, read the directory and load the betas from the underlay option.

#!/bin/bash

set -e

subjects="sub-FG01"

#subjects="sub-FG01 sub-FG02 sub-FG03 sub-FG04 sub-FG05 sub-FG06 sub-FG07 sub-FG09 sub-FG10 sub-FG11 sub-FG12 sub-FG13 sub-FG14 sub-FG16 sub-FG17 sub-FG20 sub-FG22 sub-FG23 sub-FG25 sub-FG28R sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-RG30R sub-RG31R sub-RG32R sub-RG26 sub-RG27 sub-RG30 sub-RG31 sub-RG34 sub-RG37 sub-RG38 sub-RG41 sub-RG44 sub-RG47 sub-RG48 sub-RG49 sub-RG50 sub-RG51 sub-RG55R sub-RG56R sub-RG58R sub-RG59R"

#subjects="sub-FG12 sub-FG25 sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-RG30 sub-RG31 sub-RG61R"
#working? sub-FG12 sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-RG30 sub-FG25 sub-RG61R


# select/deselect the learning phase of interest. PAV, AVO, EXT
task="Extinction"

# resample anat_roi to same resolution as master (your functional images [that you can check it via MANGO, open the func AND ROI, ctrl+I--> image dimension: are they the same??])
3dresample -master $input \
  -prefix "derivatives/ROIs/VTA_resam.nii" \
  -inset "derivatives/ROIs/VTA.nii" \
  -rmode NN
  
for subj in $subjects; do
	prefix="derivatives/afni/$subj"
	input="$prefix/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz"
	mask="derivatives/ROIs/VTA_roi_resam.nii"
	all_confounds="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_desc-confounds_timeseries.tsv"
	event_file="$subj/func/${subj}_task-${task}_events.tsv"
	averaged_BOLD_from_ROI="$prefix/averaged_BOLD_from_ROI.1D"
	
# Creates the events files (not pooled as it is for fmriprep)	
Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
X <- X[order(X$onset),] # make sure events are ordered by onset time
X <- split(X, X$trial_type) # split by trial_type
prefix <- args[2]
cat(X[["relief csm"]]$onset, sep=" ", file=file(paste0(prefix, "/relief_csm.1D")))
cat(X[["relief csav"]]$onset, sep=" ", file=file(paste0(prefix, "/relief_csav.1D")))
cat(X[["relief rating"]]$onset, sep=" ", file=file(paste0(prefix, "/relief_rating.1D")))
' $event_file $prefix

# Creates the confounds	
Rscript -e '
args = commandArgs(trailingOnly=TRUE)
X <- read.table(args[1], sep="\t", header=T, na.strings="n/a")
n <- names(X)

#CONFOUNDS 
S <- startsWith(n, "trans") | startsWith(n, "rot") | startsWith(n, "white_matter") | startsWith(n, "csf") | startsWith(n, 
X <- lapply(X[S], function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })


write.table(X, file=args[2], row.names=F, col.names=F, sep=" ")
' $all_confounds $prefix/confounds

# average voxel signal to get the mean betas using ROIs before using the events
3dmaskave -quiet -mask $mask $input > $averaged_BOLD_from_ROI

	3dDeconvolve -input $input \
	    -mask $mask \
	    -nobout \
	    -polort A \
	    -local_times \
	    -num_stimts 3 \
	    -stim_times_IM 1 "$prefix/relief_csav.1D" 'GAM' \
	    -stim_label 1 csav \
	    -stim_times_IM 2 "$prefix/relief_csm.1D" 'GAM' \
	    -stim_label 2 csm \
	    -stim_times 3 "$prefix/relief_rating.1D" 'GAM' \
	    -stim_label 3 relief_rating \
	    -ortvec $prefix/confounds "confounds" \
	    -jobs 8 \
	    -x1D $prefix/X.xmat.1D \
	    -xjpeg $prefix/X.jpg \
	    -bucket $prefix/betas_ROI \
	    -x1D_stop
	  mv 3dDeconvolve.err "$prefix/3dDeconvolve_ROI_VTA.err"
	  
	  3dREMLfit -matrix $prefix/X.xmat.1D \
	     -input ${averaged_BOLD_from_ROI}'[0]'\' \
             -Rbuck "$prefix/betas_ROI_VTA_REML" \
             -Grid 3 -tout
             
                          -Rvar "$prefix/betas_ROI_VTA_REMLvar" \
done




