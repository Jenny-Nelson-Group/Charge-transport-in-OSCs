#!/bin/bash

file=$1

for i in {0..350..10}
do

awk -f ../jkp_extract_geom.awk "${file}"_"${i}".log > "${file}"_"${i}".xyz

done
