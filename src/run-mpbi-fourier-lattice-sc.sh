#!/bin/bash
# Adjusted for Supercloud to save to #TMPDIR, which is a temporary directory on the local
# node rather than the global file system. This can speed up read/write significantly.

module load mpi/openmpi-4.0
calcname=$1
# nprocs=${2:-10}
logname=${calcname}

# odd "bug" on plasmon-ms server makes us need this: suspect this is similar to
# sec. 3.4.0.2 of https://www.quantum-espresso.org/Doc/user_guide/node19.html
#export OMP_NUM_THREADS=1
# openBLAS by default will use all available threads for BLAS operations:
# this is not helpful when we have very many cores (actually harmful)
export OPENBLAS_NUM_THREADS=1

# command substitution usually splits at spaces, which is very annoying
# here, since we need spaces in list's, vector3's etc. We disable this
# through the IFS variable (defines where command substitution will
# perform word-splitting; usually \n,\t, and space) and then reenable
# it afterwards. See https://unix.stackexchange.com/a/39482.
IFS=$'\n';
# run mpbi-mpi with inputs command-substituted from input/{$calcname}.sh
# possible extra-opts for mpirun: --report-bindings --map-by core --bind-to core
mpirun  \
    mpbi-mpi $(cat ${TMPDIR}/${calcname}.sh) \
    ctl/fourier-lattice.ctl 2>&1 | tee ${TMPDIR}/${logname}.log
# restore the value of IFS
unset IFS;

# variables
runtype=$(grep "run-type=" ${TMPDIR}/${calcname}.sh | sed 's/run-type=//;s/\"//g') # get polarization-string
res=$(grep "res=" ${TMPDIR}/${calcname}.sh | sed 's/res=//') # get resolution

# process and tidy up results
# . src/fix-unitcell.sh $calcname $runtype $res
cat ${TMPDIR}/${logname}.log | . src/get-freqs.sh $runtype ${logname}-dispersion.out ${TMPDIR}
cat ${TMPDIR}/${logname}.log | . src/get-symeigs.sh ${logname}-symeigs.out ${TMPDIR}
# . src/get-berry-phases.sh ${logname}
