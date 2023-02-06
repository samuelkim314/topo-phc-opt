#!/bin/bash

_data=$(</dev/stdin)
_pol=${1:-""}
_savename=$2
_dir=${3:-"output"}

# if _pol is "all", we actually need an empty polarization-prefix
if [ "$_pol" = "all" ]; then _pol=""; fi

echo "... saving extracted frequency data to ${_dir}/${_savename}"

echo "$_data" |
    grep "${_pol}freqs:, " |
    sed "s/${_pol}freqs:, //g; /k index.*/d; s/ //g" > \
    ${_dir}/${_savename}
