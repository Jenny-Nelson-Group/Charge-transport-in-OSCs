#!/bin/bash

file=$1

for i in {0..350..10}
do
echo $i

cat "${file}" | sed s/DIH/${i}/ > "${file%.*}"_${i}.com 

done
