#!/bin/bash

for dir in C60 M B T Tetra Penta Hexa
do

for i in {1..20}
do

OLD=" 0 50"
NEW=" $((($i*50)-50)) $(($i*50))"

sed "s/$OLD/$NEW/g" ITIAM_SCRIPT_"${dir}"_"${i}".sh > ITIAM_final_"${dir}"_"${i}".sh

done
done
