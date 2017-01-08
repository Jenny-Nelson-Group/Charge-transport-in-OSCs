#!/bin/bash

for dir in C60 #M B T Tetra Penta Hexa
do

cell=` head "${dir}"_SimulAnneal_100.pdb | grep CRYST1 | awk '{print $2,$3,$4}'`

for Time in 0 #{0..100}
do


r_vs_root_t=$(python Propagate_polaron.py Evals_unperturbedH_"${dir}"_SimulAnneal_100.pdb.xyz.dat Evecs_unperturbedH_"${dir}"_SimulAnneal_100.pdb.xyz.dat "${dir}"_SimulAnneal_100.pdb.xyz 1000 $Time "${dir}" ${cell})

echo $r_vs_root_t #>> r_vs_root_t_"${dir}"_1fs_bulk.dat

done
done

