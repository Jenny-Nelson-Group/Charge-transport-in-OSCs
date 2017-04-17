#!/bin/bash

for i in {1..50}
do
if [ -f "${i}".res ]
then

a=$(grep "CELL" "${i}".res | awk '{print $3}')  
b=$(grep "CELL" "${i}".res | awk '{print $4}')
c=$(grep "CELL" "${i}".res | awk '{print $5}')

alpha=$(grep "CELL" "${i}".res | awk '{print $6}')
beta=$(grep "CELL" "${i}".res | awk '{print $7}')
gamma=$(grep "CELL" "${i}".res | awk '{print $8}')
s=$(grep "SpaceGroup" "${i}".res)
SG=$(echo "${s##*_}" | awk '{print $1}')

echo $i $a $b $c $alpha $beta $gamma $SG >> Cell_info.dat

#sed -n -e '/C1/,$p;$d' "${i}".res | awk '{print $3,$4,$5}' > "${i}"_frac.xyz

#python Frac_to_real.py "${i}"_frac.xyz $a $b $c $alpha $beta $gamma "${i}".xyz

fi
done

