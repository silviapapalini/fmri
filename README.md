# fmri

## Pre-processing

- convert dicom to niftii+json (dcm2niigui)
- event calculation (matlab script): onset + duration of the CSs
- [BIDS validator](http://bids-standard.github.io/bids-validator/)
- fmriprep

## AFNI pre-processing

The `smooth_scale.sh` script performs some additional processing on the `fmriprep` output:

- read input from MNI152NLin2009cAsym space
- smooth 4mm
- `3dcalc` (c * min(200, a/b*100)*step(a)*step(b))
  scaled images, in which each voxel has a mean signal intensity of 100. This allows us to specify any changes relative to the mean as percent signal change; i.e., a value of 101 could be interpreted as a signal change of 1%.
  
## Pooled processing

### Individual-level

`RunDeconv_WB_ROIs_pooled_data(AVO|EXT).sh` does individual level processing for all individuals, including both whole-brain and ROI.

Steps:

- resample ROIs to dimensions of data sets
- generate stimulus files (AFNI wants separate files for each stimulus and only onset times)
- select and generate confounds

#### Whole-Brain

- `3dDeconvolve`: perform GLM given the stimuli and confounds. Produces the regressor matrix and outputs betas/tstats/fstats (`statsWB` file)
- Coefficients can be extracted from `statsWB` given the index. Use `3dinfo` to figure out which index is what.

### ROIs

Extract the ROI (averaged) betas from the WB stats. This is the same as averaging the BOLD signal first and then performing the GLM on the ROI signals.

- use `3dmaskave` with the (binary) ROI mask to average a coefficient from `statsWB`. Using its index as described above. Multiple indices can be given, each is averaged separately.

### Population-level

#### Whole-Brain

Run T-tests on the population level using `3dttest++`

- use `uber-ttest.py` to define masks,inputs and contrast
- generate script from the GUI and run (change the output folders)

#### ROIs

Run t-tests on the extracted ROI coefficients. This is performed in `R`. Not sure AFNI's `3dttest++` can use 1D data. 

First coefficients for the different contrasts need to be extracted for the ROIs and all subjects. This is pretty cumbersome because AFNI forgets the labels of the coeffients in the 1D file. But they are the same labels/order as given to the `3dmaskave` command. The `1d_to_Csv` script loops over all subjects and ROIs to collect the coeffients and relabel them correctly.

This produced CSV file can then be read into `R` and t-test can be done.
