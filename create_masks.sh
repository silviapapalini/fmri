#!/bin/bash

set -e

recompute=true
subjects="sub-FG01 sub-FG02 sub-FG03 sub-FG04 sub-FG05 sub-FG06 sub-FG07 sub-FG09 sub-FG10 sub-FG11 sub-FG12 sub-FG14 sub-FG15 sub-FG16 sub-FG17 sub-FG20 sub-FG22 sub-FG23 sub-FG25 sub-FG28R sub-FG29R sub-FG30R sub-FG31R sub-FG32R sub-FG33R sub-RG26 sub-RG27 sub-RG30 sub-RG31 sub-RG32 sub-RG33 sub-RG34 sub-RG35 sub-RG37 sub-RG38 sub-RG41 sub-RG44 sub-RG47 sub-RG48 sub-RG49 sub-RG50 sub-RG51 sub-RG55R sub-RG56R sub-RG58R sub-RG59R sub-RG61R sub-RG62R sub-RG63R sub-RG64R"

VTA_mask=derivatives/ROIs/canlab_VTA_bilateral.nii
 
smooth_input=$(mktemp -u --suffix=.nii)
function smooth_and_threshold_mask {
   input=$1
   output=$2
   
   [ -e $output -a $recompute == "true" ] && rm $output
   if [ ! -e $output ]; then
	   3dmerge -1blur_fwhm 4.0 -doall -prefix $smooth_input $input && \
	   3dcalc -a ${smooth_input} -expr "step(a-0.15)" -prefix ${output} &&
	   rm $smooth_input
   fi
}

function threshold_mask {
   input=$1
   output=$2
   
   [ -e $output -a $recompute == "true" ] && rm $output
   if [ ! -e $output ]; then
	   3dcalc -a ${input} -expr "step(a-0.15)" -prefix ${output}
   fi
}

function include_vta {
   input=$1
   output=$2
   
   [ -e $output -a $recompute == "true" ] && rm $output
   if [ ! -e $output ]; then
           tmpfile=$(mktemp -u --suffix=.nii)
	   3dresample -master ${input} -prefix ${tmpfile} -inset ${VTA_mask} -rmode NN   
	   3dmask_tool -input ${input} ${tmpfile} -prefix ${output} -union
	   rm $tmpfile
   fi
}

for subj in $subjects; do
	prefix="derivatives/afni/$subj"
	gm_mask="${prefix}/gm_mask.nii.gz"
	gm_vta_mask="${prefix}/gm_vta_mask.nii.gz"
	mask_dir="derivatives/fmriprep/$subj/anat"
	gm_prob="${mask_dir}/${subj}_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz"
	
	smooth_and_threshold_mask $gm_prob $gm_mask
	include_vta $gm_mask $gm_vta_mask
done

function all_masks {
    mask_name=$1
    shift

    array=()
    for arg in "$@"; do
        array+=("derivatives/afni/$arg/$mask_name")
    done

    echo "${array[@]}"
}

masks=$(all_masks gm_mask.nii.gz $subjects)

mean_probseg="derivatives/afni/SL/mean_gm_probseg.nii.gz"
[ -e $mean_probseg -a $recompute == "true" ] && rm $mean_probseg
[ ! -e ${mean_probseg} ] && 3dMean -prefix ${mean_probseg} $masks

mean_mask="derivatives/afni/SL/mean_gm_mask.nii.gz"
[ -e $mean_mask -a $recompute == "true" ] && rm $mean_mask
3dcalc -a ${mean_probseg} -expr 'step(a-0.04)' -prefix $mean_mask

mean_vta_mask="derivatives/afni/SL/mean_gm_vta_mask.nii.gz"
include_vta $mean_mask $mean_vta_mask

