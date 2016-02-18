#!/bin/bash

for dir in C60 M B T Tetra Penta Hexa
do

#mv Evals_localising_all_sites_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_0.dat Evals_localising_all_sites_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_000.dat

#mv Evecs_localising_all_sites_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_0.dat Evecs_localising_all_sites_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_000.dat

#mv Perturbed_H_diag_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_0.dat Perturbed_H_diag_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_000.dat

#mv Evals_localising_all_sites_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_50.dat Evals_localising_all_sites_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_050.dat

#mv Evecs_localising_all_sites_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_50.dat Evecs_localising_all_sites_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_050.dat

#mv Perturbed_H_diag_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_50.dat Perturbed_H_diag_"${dir}"_SimulAnneal_100.pdb.xyz_0.32_050.dat

cat Evecs_localising_all_sites_"${dir}"_* > Evecs_"${dir}"_all.dat

cat Evals_localising_all_sites_"${dir}"_* > Evals_"${dir}"_all.dat

cat Perturbed_H_diag_"${dir}"_* > Perturbed_H_diag_"${dir}"_all.dat

done 
