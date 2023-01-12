#!/bin/bash
_calcname=$1
_runtype=$2
_res=$3

epsfilebase=$(echo $_calcname | sed "s/-${_runtype}//")

# check if a filename (or multiple!) $_calcname*-epsilon.h5 exists; if so, process 

# unit cell and save result to output/unitcell/${epsfilebase}-epsilon.h5
for foundh5epsfile in ${_calcname}*-epsilon.h5; do # wildcard if-exists-then
    [ -f "$foundh5epsfile" ] || continue # quit iteration if foundh5epsfile is not a file 
    echo ""; echo "... expanding epsilon .h5 file to cartesian basis, saving to:"
    mpb-data -r -m 1 -n $_res $foundh5epsfile
    mv -f -v $foundh5epsfile "output/unitcell/${epsfilebase}-epsilon.h5"
done
