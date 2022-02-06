smooth_output=$(mktemp -u --suffix=.nii)
mean_output=$(mktemp -u --suffix=.nii)

for subj in `cat participants.tsv | tr -d '\r'` ; do
	if [ -d "derivatives/fmriprep/$subj/func" ]; then
	  echo "generating for $subj"
	  mkdir -p "derivatives/afni/${subj}"
	  for task in Avoidance Extinction Pavlovian; do
	     input="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"
	     mask="derivatives/fmriprep/${subj}/func/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"
	     output="derivatives/afni/${subj}/${task}/${subj}_task-${task}_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_scaled.nii.gz"
	     if [ ! -e $output -o $input -nt $output ]; then
	       [ -e $output ] && rm $output
	       echo "creating $output"
	       3dmerge -1blur_fwhm 4.0 -doall -prefix $smooth_output $input &&
	       3dTstat -prefix $mean_output $smooth_output &&
	       3dcalc -a $smooth_output -b $mean_output -c $mask \
	       	-expr 'c * min(200, a/b*100)*step(a)*step(b)' \
	       	-prefix $output &&
	       rm $smooth_output $mean_output
	     fi
	  done
	fi
done


