#!/bin/bash

function label2index {
   if [ "$#" -ne 2 ]; then
      echo "label2index: expects dataset and label" >&2
      return 1
   fi
   
   dset=$1
   label=$2
   index=$(3dinfo -label2index $label $dset) || return 2

   if [ -z "$index" ]; then
      echo "label2index: failed to get index for $label" >&2
      return 3
   fi

   echo $index
}

# This function creates a sphere around a coordinate xyz (LPI)
function make_spherical_roi {
   if [ "$#" -ne 4 ]; then
      echo "label2index: <master> <xyz> <prefix> <radius>" >&2
      return 1
   fi
   
   master=$1
   xyz=$2
   prefix=$3
   radius=$4

   echo $xyz | 3dUndump -orient LPI -srad $radius -master $master -prefix $prefix -xyz -
}

# TESTS
#index=$(label2index $1 $2) && echo "Index found is $index"
