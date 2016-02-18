#!/bin/bash

for dir in C60 M B T Tetra Penta Hexa
do
for i in {1..20}
do
qsub ITIAM_script_"${dir}"_"${i}".sh
done
done
