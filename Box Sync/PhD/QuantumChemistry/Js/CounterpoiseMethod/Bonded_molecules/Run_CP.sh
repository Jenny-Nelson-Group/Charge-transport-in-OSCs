#!/bin/bash

file=$1

deghomo=$2
deglumo=$3

for i in {0..350..10}
do

echo $i

python ../Counterpoise_method.py "${file}"_"${i}"A.log "${file}"_"${i}"B.log "${file}"_"${i}"AB.log "${deghomo}" "${deglumo}" > J_conv_"${file}"_"${i}".txt

grep "HOMO-HOMO coupling" J_"${file}"_"${i}".txt > temph.txt
grep "LUMO-LUMO coupling" J_"${file}"_"${i}".txt > templ.txt

J_homo=$(grep -Eo '\<[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\>' temph.txt)
J_lumo=$(grep -Eo '\<[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\>' templ.txt)
 
echo $i $J_homo $J_lumo >> J_change_dihedral_"${file}".txt

done



