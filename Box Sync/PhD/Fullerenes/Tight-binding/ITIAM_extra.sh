#!/bin/sh 

pdb=$1 
SCstart=$2
SCend=$3


for pdb in $1
do
 cell=` head "${pdb}" | grep CRYST1 | awk '{print $2,$3,$4}' ` # extract periodic cell vectors
 grep "^ATOM" "${pdb}" | grep " C " | awk '{print $6,$7,$8}' > "${pdb}".xyz #all the 'C'60 pseudo-atoms
 sites=` cat "${pdb}".xyz | wc -l ` #cat'ing so just have single field of lines


for alpha in 0.2 0.32 #0.001 0.01 0.05 0.1 0.32
do
for dx in 0 #0.001 0.005 0.01 0.05 0.1
do
for state in 0
do

if [ ! -f $WORK/Fullerenes/Polarons/Evals_"${pdb}".xyz_"${alpha}"_"${dx}"_"${SCstart}".dat ]

then
echo ${cell}
echo $alpha

python ITIAM_final.py ${sites} "${pdb}".xyz ${cell} $alpha $dx $state $SCstart $SCend

fi
done
done
done



done



