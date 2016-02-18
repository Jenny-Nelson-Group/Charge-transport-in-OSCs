#!/bin/bash

for dir in C60 M B T Tetra Penta Hexa
do
for num in 0{91..93}
do
for i in {1..20}
do
qsub ITIAM_final_"${dir}"_"${num}"_"${i}".sh
done
done
done
