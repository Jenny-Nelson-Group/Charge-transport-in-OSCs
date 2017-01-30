#!/bin/bash

file=$1
namefile=$2

for i in 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0 #{0..350..10}
do

bash ../Bonded_counterpoise.sh "${file}"_"${i}".xyz "${namefile}"_"${i}"

done

