#!/usr/bin/bash

cd D
../sort.sh D
../sum.py D
../dens_rel_all.plt D 2
cd ../He
../sort.sh He
../sum.py He
./dens_rel_all.plt He 4
cd ../Be
../sort.sh Be
../sum.py Be
../dens_rel_all.plt Be 9
cd ../C
../sort.sh C
../sum.py C
./dens_rel_all.plt C 12
cd ../O
../sort.sh O
../sum.py O
../dens_rel_all.plt O 16

cd ..
okular */*.pdf &
#cd ../Al
#./sum.py Al
#./dens_rel_all.plt Al 27
#cd ../Fe
#./sum.py Fe
#./dens_rel_all.plt Fe 56

