#!/bin/bash

for i in 0{1..9} {10..49}
do
if [ -f "${i}".cif ]
then

babel -icif "${i}".cif -oxyz "${i}".xyz

fi
done

