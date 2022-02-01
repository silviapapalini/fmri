# fmri

## Pre-processing

- convert dicom to niftii+json (dcm2niigui)
- event calculation (matlab script): onset + duration of the CSs
- [BIDS validator](http://bids-standard.github.io/bids-validator/)
- fmriprep

## AFNI

- input from MNI152NLin2009cAsym space
- smooth 4mm
- 3dcalc (c * min(200, a/b*100)*step(a)*step(b))
  scaled images, in which each voxel has a mean signal intensity of 100. This allows us to specify any changes relative to the mean as percent signal change; i.e., a value of 101 could be interpreted as a signal change of 1%.
