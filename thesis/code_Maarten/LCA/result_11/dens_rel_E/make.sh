#!/usr/bin/bash

cd D
./dens_rel_all.plt D 2
cd ../He
./dens_rel_all.plt He 4
cd ../C
./dens_rel_all.plt C 12
cd ../Al
./sum.py Al
./dens_rel_all.plt Al 27
cd ../Fe
./sum.py Fe
./dens_rel_all.plt Fe 56

