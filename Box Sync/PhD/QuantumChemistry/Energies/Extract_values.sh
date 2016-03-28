#!/bin/bash

file=$1
Type=$2

for i in 30_31 31_30 20_38 24_25 52_35 7_8  #10_26 11_12 13_14 15_32 16_17 18_19 20_38 21_22 24_25 27_28 29_47 30_31 33_34 35_52 36_37 39_40 41_42 43_44 45_46 48_49 50_51 53_54 55_56 5_6 57_58 59_60 61_62 7_8 9_23
do

if [ "$Type" == "HOMO" ];
then
 grep "occ" "${file}"_"${i}"_neutral.log | tail -n 1 | awk '{print ($NF*27.211)}'
fi

if [ "$Type" == "LUMO" ];
then
 grep "virt" "${file}"_"${i}"_neutral.log | head -n 1 | awk '{print $5*27.211 }'
fi

if [ "$Type" == "SCF_HOMO" ];
then
neutral=` grep "SCF Done" "${file}"_"${i}"_neutral.log | awk '{print $5}' `
cation=` grep "SCF Done" "${file}"_"${i}"_cation.log | awk '{print $5}' `
echo "($neutral - $cation ) * 27.211" | bc -l
fi

if [ "$Type" == "SCF_LUMO" ];
then
neutral=` grep "SCF Done" "${file}"_"${i}"_neutral.log | awk '{print $5}' `
anion=` grep "SCF Done" "${file}"_"${i}"_anion.log | awk '{print $5}' `
echo "($anion - $neutral ) * 27.211" | bc -l
fi

if [ "$Type" == "Energy" ];
then
neutral=` grep "SCF Done" "${file}"_"${i}"_neutral.log | awk '{print $5}' `
echo $neutral
fi


done


