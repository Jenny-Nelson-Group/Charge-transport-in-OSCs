#!/bin/bash

for i in 0{1..9} {10..49}
do
if [ -f "${i}".cif ]
then

x=$(grep "_cell_length_a" "${i}".cif | awk '{print $2}')  
y=$(grep "_cell_length_b" "${i}".cif | awk '{print $2}')
z=$(grep "_cell_length_c" "${i}".cif | awk '{print $2}')

echo $i $x $y $z >> Cell_lengths.dat

fi
done

