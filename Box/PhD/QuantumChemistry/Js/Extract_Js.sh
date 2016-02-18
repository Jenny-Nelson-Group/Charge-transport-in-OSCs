#!/bin/sh

dir=$1
mol1=$2
mol2=$3


for i in "${mol1}"  
do
for j in "${mol2}"
do
if [ -f ""${dir}"_"$i"_"$j"_J_nosymm.txt" ];
then

grep "Non_degenerate_LUMO_coupling" "${dir}"_"$i"_"$j"_J_nosymm.txt > temp.txt

J=$(grep -Eo '\<[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\>' temp.txt) 

echo $i $j $J >> Js_LUMO_"${dir}".txt

#else 

#echo $i $j 0 >> Js_HOMO_EP.txt 

fi

done
done
