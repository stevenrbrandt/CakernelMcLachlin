#! /bin/bash

MATHEMATICA="/Applications/Mathematica.app/Contents/MacOS/MathKernel"

# Run Mathematica to regenerate the code
< McLachlan.m "$MATHEMATICA" | tee McLachlan.out
