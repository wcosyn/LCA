#!/usr/bin/bash

cd H
./plot_ob.plt
./plot_xx_ob.plt
cd ../He
./plot_ob.plt
./plot_xx_ob.plt
../plot_xx_ratio.plt He 4
cd ../Be
./plot_ob.plt
#./plot_xx_ob.plt
#../plot_xx_ratio.plt Be 9
cd ../C
./plot_ob.plt
./plot_xx_ob.plt
../plot_xx_ratio.plt C 12
cd ../O
./plot_ob.plt
#./plot_xx_ob.plt
#../plot_xx_ratio.plt O 16
cd ../Al
./plot_ob.plt
./plot_xx_ob.plt
../plot_xx_ratio.plt Al 27
cd ../Ca40
./plot_ob.plt
./plot_xx_ob.plt
../plot_xx_ratio.plt Ca 40
cd ../Ca48
./plot_ob.plt
./plot_xx_ob.plt
../plot_xx_ratio.plt Ca 48
cd ../Fe
./plot_ob.plt
./plot_xx_ob.plt
../plot_xx_ratio.plt Fe 56
cd ../Ag
./plot_ob.plt
#./plot_xx_ob.plt
#../plot_xx_ratio.plt Ag 108
