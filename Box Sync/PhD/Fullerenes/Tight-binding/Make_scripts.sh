#!/bin/bash

for dir in C60 M B T Tetra Penta Hexa
do
for i in {1..20}
do

OLD="ITIAM_everysite.sh"
NEW="ITIAM_final.sh"

#cp ITIAM_script_C60_"${i}".sh ITIAM_script_"${dir}"_"${i}".sh

sed "s/$OLD/$NEW/g" ITIAM_script_"${dir}"_"${i}".sh > ITIAM_SCRIPT_"${dir}"_"${i}".sh  

done
done
