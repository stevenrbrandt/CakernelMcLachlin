#! /bin/bash

name=$1

# Create a thorn $name_Helper

orig=ML_BSSN

cp -R prototype/${orig}_Helper ${name}_Helper
find ${name}_Helper -name '*~' | xargs rm
find ${name}_Helper -type f | xargs perl -pi -e "s/${orig}/${name}/g"
