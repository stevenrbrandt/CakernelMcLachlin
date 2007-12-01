#! /bin/bash

# Abort on errors
set -e

MATHEMATICA="math"

script=$1

if test -z "$script"; then
    echo "Usage:"
    echo "$0 <script.m>"
    exit 2
fi

error=$(basename $script .m).err
output=$(basename $script .m).out

rm -f $output

# Run Mathematica to regenerate the code
< $script "$MATHEMATICA" | tee $error

if grep 'KrancError' $error; then
    echo
    echo "There was an error when running Kranc on $script."
    echo "The file $error contains details."
    echo
    echo "*** The Cactus thorns have NOT been updated. ***"
    echo
    exit 1
fi

mv $error $output
