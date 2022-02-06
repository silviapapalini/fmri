#!/bin/bash
afni="derivatives/afni/"
if [ ! -d $afni ]; then
  echo "AFNI folder not found: $afni" >&2
  exit 1;
fi

ls $afni/sub-*/3dDeconvolve.err $afni/sub-*/stats* $afni/sub-*/errts* $afni/sub-*/fitts* $afni/sub-*/X.* $afni/sub-*/confounds

read -p "Are you sure to delete the files above? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]; then
   rm $afni/sub-*/3dDeconvolve.err $afni/sub-*/stats* $afni/sub-*/errts* $afni/sub-*/fitts* $afni/sub-*/X.* $afni/sub-*/confounds
fi
