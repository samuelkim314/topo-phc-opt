#!/bin/bash

calcname=$1

# openBLAS by default will use all available threads for BLAS operations:
# this is not helpful when we have very many cores (actually harmful)
export OPENBLAS_NUM_THREADS=1 

# command substitution usually splits at spaces, which is very annoying
# here, since we need spaces in list's, vector3's etc. We disable this
# through the IFS variable (defines where command substitution will 
# perform word-splitting; usually \n,\t, and space) and then reenable 
# it afterwards. See https://unix.stackexchange.com/a/39482.
IFS=$'\n'; 
# run mpb with inputs command-substituted from input/{$calcname}.sh
~/postdoc/mpb-transform-dev/1.8-dev/bin/mpb \
    $(cat input/${calcname}.sh) \
    berry-bool="true" \
    ctl/fourier-lattice.ctl 2>&1 | tee logs/${calcname}-berry.log
# restore the value of IFS
unset IFS;

# variables
runtype=$(grep "run-type=" input/${calcname}.sh | sed 's/run-type=//;s/\"//g') # get polarization-string
res=$(grep "res=" input/${calcname}.sh | sed 's/res=//') # get resolution

# process and tidy up results
. fix-unitcell.sh $calcname $runtype $res
cat logs/${calcname}-berry.log | 
    sed -n '/berry phases:, /q;p' | # filter on frequencies calculated _before_ berry-phase calculation
    get-freqs.sh $runtype ${calcname}-berry-dispersion.out
cat logs/${calcname}-berry.log | 
    sed -n '/berry phases:, /,$p' | # filter on frequencies calculated _after_ berry-phase calculation
    get-freqs.sh $runtype ${calcname}-berry-loopfreqs.out
. get-berry-phases.sh ${calcname}-berry
