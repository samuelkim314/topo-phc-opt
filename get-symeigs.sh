#!/bin/bash
_data=$(</dev/stdin)
_savename=$1
_dir=${2:-"output"}

_symeigs_data=$(echo "$_data" | grep "sym-eigs:, " | sed -E 's/\(|\)/"/g; s/sym-eigs:, //g; s/ //g')

# no need to save anything if there are no "sym-eigs:, " matches
if test ! -z "$_symeigs_data" # test whether _symeigs_data is a string of nonzero length
then 
    echo "... saving extracted symmetry eigenvalues data to ${_dir}/$_savename";
    echo "$_symeigs_data" > ${_dir}/$_savename
else
    echo "... no associated symmetry data to save"
fi
