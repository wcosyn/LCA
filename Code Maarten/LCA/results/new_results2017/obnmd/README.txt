

the one-body momentum distribution results
===========================================

naming convention:

dens_ob4.[xy].[uvw].[A].[abcd](s)

x : isospin selection of particle "1" in pair (1: proton, -1: neutron, 0: no sel.)
y : isospin selection of particle "2" in pair (1: proton, -1: neutron, 0: no sel.)
u : central  correlation (0 inactive, 1 active)
v : tensor   correlation (0 inactive, 1 active)
w : spin/iso correlation (0 inactive, 1 active)
a : n quantum number selection of rcm. state "1" (-1: no selection)
b : l quantum number selection of rcm. state "1" (-1: no selection)
c : n quantum number selection of rcm. state "2" (-1: no selection)
d : l quantum number selection of rcm. state "2" (-1: no selection)
s : (optional) sorted file. Because multithreading files are generally not
sorted in one-body momentum k. (sort script present in this dir)
