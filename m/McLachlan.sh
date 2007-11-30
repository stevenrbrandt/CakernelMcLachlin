#! /bin/bash

# Run this script in the "m" subdirectory.
# It re-generates the Cactus thorns trees and copies them if they have
# changed.

# Abort on errors
set -e

MATHEMATICA="math"

THORNS="ML_ADM ML_BSSN ML_WaveToy"
SCRIPTS="McLachlan.m WaveToy.m"

# Remove old output
rm -rf $THORNS

for script in $SCRIPTS; do
    output=$(basename $script .m).out
    
    # Run Mathematica to regenerate the code
    < $script "$MATHEMATICA" | tee $output
    
    if grep 'KrancError' $output > /dev/null 2>&1; then
        echo
        echo "There was an error when running Kranc on $script."
        echo "The file $output contains details."
        echo
        echo "*** The Cactus thorns have NOT been updated. ***"
        echo
        exit 1
    fi
done

# Copy the source trees
for thorn in $THORNS; do
    ./copy-if-changed.sh $thorn ../$thorn
done

echo
echo "The Cactus thorns have been regenerated successfully."
echo
