#!/bin/bash

file=$1
dir=$2

for Type in C60 #PCBM
do
for field in 100 #316 3162 1000 9999
do


grep "${Type}"_.p0.log "${file}" > "${Type}"_p0.txt
grep "${Type}"_.p"${field}".log "${file}" > "${Type}"_p"${field}".txt

python ../Calculate_epsilon.py "${Type}"_"${dir}"p"${field}".txt "${Type}"_"${dir}"p0.txt
   
done
done

