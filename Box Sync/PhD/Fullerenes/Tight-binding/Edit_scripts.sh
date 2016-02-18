#!/bin/bash

for dir in C60 M B T Tetra Penta Hexa
do
for i in {1..20}
do

OLD="walltime=3:58:02"
NEW="walltime=20:58:02"


sed "s/$OLD/$NEW/g" ITIAM_final_"${dir}"_"${i}".sh > ITIAM_final_"${dir}"_091_"${i}".sh  

done
done
