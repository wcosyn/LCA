#!/usr/bin/bash

cd D
./dens_rel_pair H 2
cd ../He
./dens_rel_pair.plt He 4
cd ../C
./dens_rel_pair.plt C 12
cd ../Al
./dens_rel_pair.plt Al 27
cd ../Fe
./dens_rel_pair.plt Fe 56

