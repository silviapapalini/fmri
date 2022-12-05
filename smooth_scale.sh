set -e

wm_inv_mask=derivatives/ROIs/WM_inv_mask.nii.gz
ROIs="VTA canlab_VTA_bilateral NAcc canlab_NAC_bilateral VmPFC VmPFC_box VmPFC_csp_vs_csm"

mean_output=$(mktemp -u --suffix=.nii)
roi_mean=$(mktemp -u --suffix=.1D)
intermediate_tmp=$(mktemp -u --suffix=.nii)

for subj in `cat participants.tsv | tr -d '\r'` ; do
	if [ -d "derivatives/fmriprep/$subj/func" ]; then
	  echo "generating for $subj"
	  mkdir -p "derivatives/afni/${subj}"
	  prefix="derivatives/afni/$subj"
	  
	  master="derivatives/fmriprep/${subj}/func/${subj}_task-Avoidance_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"
	  	  
  	  # resample WM_inv
	  if [ ! -e "$prefix/WM_inv_resam.nii" ]; then
		3dresample -master $master \
		  -prefix "$prefix/WM_inv_resam.nii" \
		  -inset $wm_inv_mask \
		  -rmode NN
	  fi

	  # Ceeate ROIs for subject
	  if [ ! -e "$prefix/VTA_resam.nii" ]; then
		3dresample -master $master \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/VTA_bram.nii.gz" \
		  -rmode NN
		mv $intermediate_tmp "$prefix/VTA_resam.nii" # original VTA from Esser (VTA plus SN)
	  fi
	  if [ ! -e "$prefix/canlab_VTA_bilateral_resam.nii" ]; then
		3dresample -master $master \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/canlab_VTA_bilateral.nii" \
		  -rmode NN
		mv $intermediate_tmp "$prefix/canlab_VTA_bilateral_resam.nii" # original VTA from Esser (VTA plus SN)
	  fi


	  if [ ! -e "$prefix/NAcc_resam.nii" ]; then
		3dresample -master $master \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/NAcc_HarvardOxford.nii.gz" \
		  -rmode NN

		3dmask_tool -input "$prefix/WM_inv_resam.nii" $intermediate_tmp -prefix "$prefix/NAcc_resam.nii" -inter
		rm $intermediate_tmp
	  fi
	  if [ ! -e "$prefix/canlab_NAC_bilateral_resam.nii" ]; then
		3dresample -master $master \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/canlab_NAC_bilateral.nii" \
		  -rmode NN

		mv $intermediate_tmp "$prefix/canlab_NAC_bilateral_resam.nii"
	  fi  

	  if [ ! -e "$prefix/VmPFC_resam.nii" ]; then
		3dresample -master $master \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/VmPFC_parcels.nii.gz" \
		  -rmode NN

		3dmask_tool -input "$prefix/WM_inv_resam.nii" $intermediate_tmp -prefix "$prefix/VmPFC_resam.nii" -inter
		rm $intermediate_tmp
	  fi	  
	  
	  if [ ! -e "$prefix/VmPFC_box_resam.nii" ]; then
		3dresample -master $master \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/VmPFC_box.nii" \
		  -rmode NN

		mv $intermediate_tmp "$prefix/VmPFC_box_resam.nii"
	  fi		
	  
	  if [ ! -e "$prefix/VmPFC_csp_vs_csm_resam.nii" ]; then
		3dresample -master $master \
		  -prefix $intermediate_tmp \
		  -inset "derivatives/ROIs/VmPFC_csp_vs_csm.nii" \
		  -rmode NN

		mv $intermediate_tmp "$prefix/VmPFC_csp_vs_csm_resam.nii"
	  fi  
	  
	  
	  for task in Avoidance Extinction Pavlovian; do
	     input="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"
	     mask="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"
	     smooth_output="derivatives/afni/${subj}/${task}/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth.nii.gz"
	     scaled_output="derivatives/afni/${subj}/${task}/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz"
	     if [ -e $input ]; then
  	       if [ ! -e $smooth_output -o $input -nt $smooth_output ]; then
	         [ -e $smooth_output ] && rm $smooth_output
	         echo "creating $smooth_output"
	         3dmerge -1blur_fwhm 4.0 -doall -prefix $smooth_output $input
	         if [ ! -e $scaled_output ]; then #-o $smooth_output -nt $scaled_output ]; then
  	           3dTstat -prefix $mean_output $smooth_output
 	           3dcalc -a $smooth_output -b $mean_output -c $mask \
	       	-expr 'c * min(200, a/b*100)*step(a)*step(b)' \
	       	-prefix $scaled_output
	           rm $mean_output
	         fi
	       fi
	     
	       for roi in $ROIs; do
	         roi_mask="${prefix}/${roi}_resam.nii"
	         averaged_BOLD_from_ROI="derivatives/afni/${subj}/${task}/${subj}_task-${task}_averaged_BOLD_from_${roi}.1D"
	         scaled_BOLD_from_ROI="derivatives/afni/${subj}/${task}/${subj}_task-${task}_scaled_BOLD_from_${roi}.1D"
	         if [ ! -e $scaled_BOLD_from_ROI ]; then
  	           3dmaskave -quiet -mask $roi_mask $smooth_output > $averaged_BOLD_from_ROI
	           3dTstat -prefix $roi_mean $averaged_BOLD_from_ROI\'
	           3dcalc -a $averaged_BOLD_from_ROI\' -b $roi_mean -expr 'min(200, a/b*100)*step(a)*step(b)' -prefix $scaled_BOLD_from_ROI
	           rm $roi_mean
	         fi
	       done
	    fi
	  done
	fi
done


