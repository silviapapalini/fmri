# fmri

## Pre-processing

- convert dicom to niftii+json (dcm2niigui)
- event calculation (matlab script): onset + duration of the CSs
- [BIDS validator](http://bids-standard.github.io/bids-validator/)
- fmriprep

## AFNI

- blur 4mm
- 3dcalc (c * min(200, a/b*100)*step(a)*step(b))
