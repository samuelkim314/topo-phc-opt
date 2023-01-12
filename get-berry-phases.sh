#!/bin/bash

_calcname=$1

# look for berry phases data
_berry_data=$(grep "berry phases:, " logs/${_calcname}.log | sed 's/berry phases:, //g; s/ //g' | tail -n +2)
#grep "berry phases:, " logs/${_calcname}.log \
#    | sed 's/berry phases:, //g; s/ //g' | tail -n +2 \
#    | tee output/${_calcname}-phases.out

# no need to save anything if there are no "berry phases:, " matches
if test ! -z "$_berry_data" # test whether _berry_data is a string of nonzero length
then 
    echo "... saving extracted berry phases to output/${_calcname}-phases.out";
    echo "$_berry_data" > output/${_calcname}-phases.out

    echo "... saving associated berry multiplet structure to output/${_calcname}-multiplet-structure.out";
    grep "berry phases:, " logs/${_calcname}.log | head -n 1 \
        | sed 's/berry phases:, loop-start k1, loop-start k2, loop-start k3, loop G1, loop G2, loop G3, //; s/band (//g; s/multiplet (//g; s/)/\n/g; s/, //g' \
        | sed '/^$/d' \
        | tee output/${_calcname}-multiplet-structure.out
else
    echo "... no associated berry phase data to save"
fi
