#!/bin/sh

dir=$1
mol1=$2
mol2=$3

for i in "${mol1}"
do
for j in "${mol2}"
do
#if [ $i -lt $j ];
#then
#if [ -f ""${dir}"_"$i"_"$j"_J.log" ];
#then

echo "${dir}"_"$i"_"$j".log 

grep "basis functions," ""${dir}"_"$i"_J.log" | awk '{print $"1"}'>tmp
grep "alpha electrons" ""${dir}"_"$i"_J.log" | awk '{print $"1"}'>>tmp

grep "basis functions," ""${dir}"_"$j"_J.log" | awk '{print $"1"}'>>tmp
grep "alpha electrons" ""${dir}"_"$j"_J.log" | awk '{print $"1"}'>>tmp

grep "basis functions," ""${dir}"_"$i"_"$j"_J.log" | awk '{print $"1"}'>>tmp
grep "alpha electrons" ""${dir}"_"$i"_"$j"_J.log" | awk '{print $"1"}'>>tmp

echo "temp file made"

./rewrite_S_phi_E.x ""${dir}"_"$i"_J" ""${dir}"_"$j"_J" ""${dir}"_"$i"_"$j"_J"

cat tmp | octave get_J_interactive.m > ""${dir}"_"$i"_"$j"_J_nosymm.txt"  

echo "J calculated for molecules " $i "and" $j  

#fi
#fi
done
done

