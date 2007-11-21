#! /bin/bash

# Run this script in the "m" subdirectory.
# It re-generates the Cactus thorns trees and copies them if they have
# changed.

# Abort on errors
set -e

MATHEMATICA="math"

# Remove old output
rm -rf ML_ADM ML_BSSN

# Run Mathematica to regenerate the code
< McLachlan.m "$MATHEMATICA" | tee McLachlan.out

if grep 'KrancError' McLachlan.out > /dev/null 2>&1; then
    echo
    echo "There was an error when running Kranc."
    echo "The file McLachlan.out contains details."
    echo
    echo "*** The Cactus thorns have NOT been updated. ***"
    echo
    exit 1
fi

# Copy the source trees
./copy-if-changed.sh ML_ADM ../ML_ADM
./copy-if-changed.sh ML_BSSN ../ML_BSSN

echo
echo "The Cactus thorns have been regenerated successfully."
echo
