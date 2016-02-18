#!/bin/bash

for dir in C60 #M B T Tetra Penta Hexa
do
for n in 100
do
for Type in Evals Evecs Perturbed_H_diagonal
do
for alpha in 0.001 0.01 0.05 0.1 #0.32 
do
for dx in 0.0 #0.001 0.005 0.01 0.05 0.1
do

if [ -f "${Type}"_"${dir}"_SimulAnneal_"${n}".pdb.xyz_"${alpha}"_"${dx}"_0.dat ]
then

mv "${Type}"_"${dir}"_SimulAnneal_"${n}".pdb.xyz_"${alpha}"_"${dx}"_0.dat "${Type}"_"${dir}"_SimulAnneal_"${n}".pdb.xyz_"${alpha}"_"${dx}"_000.dat

fi

if [ -f "${Type}"_"${dir}"_SimulAnneal_"${n}".pdb.xyz_"${alpha}"_"${dx}"_50.dat ]
then

mv "${Type}"_"${dir}"_SimulAnneal_"${n}".pdb.xyz_"${alpha}"_"${dx}"_50.dat "${Type}"_"${dir}"_SimulAnneal_"${n}".pdb.xyz_"${alpha}"_"${dx}"_050.dat

fi

for num in 000 050 {100..950..50}
do

sed 's/ /\n/g' "${Type}"_"${dir}"_SimulAnneal_"${n}".pdb.xyz_"${alpha}"_"${dx}"_"${num}".dat > "${Type}"_"${dir}"_"${n}"_"${alpha}"_"${dx}"_"${num}"_polaron.dat

done

cat "${Type}"_"${dir}"_"${n}"_"${alpha}"_"${dx}"_*_polaron.dat > "${Type}"_"${dir}"_"${n}"_"${alpha}"_"${dx}"_polarons.dat

done
done
done
done
done
