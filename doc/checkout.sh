#! /bin/bash

set +e
set +x

# Cactus
#wget http://preview.cactuscode.org/download/GetCactus
echo ':pserver:cvs_anon@cvs.cactuscode.org:/cactus Ay=0=' > cvspass
env CVS_PASSFILE=cvspass cvs -d :pserver:cvs_anon@cvs.cactuscode.org:/cactus checkout -d qqq Utilities/Scripts/GetCactus
mv qqq/GetCactus .
rmdir -rf qqq
#cp /Users/eschnett/Utilities/Scripts/GetCactus .
chmod a+x GetCactus
{ echo; echo; echo 2; echo; } | ./GetCactus

# Basic thorns
pushd Cactus
{ echo; echo; echo q; } | gmake checkout
popd

# Carpet
pushd Cactus
git clone -o carpet git://carpetcode.dyndns.org/carpet.git
cd arrangements
ln -s ../carpet/Carpet* .
popd

# McLachlan
pushd Cactus
cd arrangements
git clone git://carpetcode.dyndns.org/McLachlan.git
popd

# Kranc
pushd Cactus
git clone http://www.aei.mpg.de/~ianhin/kranc.git
cd arrangements
ln -s ../kranc/Auxiliary/Cactus/KrancNumericalTools .
popd

# All other public thorns (AEIThorns, AstroGrid, LSUThorns, TAT,
# Whisky)
{ echo; echo; echo; echo; echo; } | ./GetCactus qc0-mclachlan-public.th
