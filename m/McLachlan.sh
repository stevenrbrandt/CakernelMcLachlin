#! /bin/bash

# Run this script in the "m" subdirectory.
# It re-generates the Cactus thorns trees and copies them if they have
# changed.

MATHEMATICA="/Applications/Mathematica.app/Contents/MacOS/MathKernel"

# Run Mathematica to regenerate the code
< McLachlan.m "$MATHEMATICA" | tee McLachlan.out

# Copy the source trees
./copy-if-changed.sh ML_ADM ../ML_ADM
./copy-if-changed.sh ML_BSSN ../ML_BSSN
