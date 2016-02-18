#!/bin/sh
#  runDoS_from_CGPCBM_simulations.sh - generates Js (.edges) file with ./cg-pcbm-j for specified PDB PCBM CG files
#       & runs python TB ToS calculator (DBTW.py)
#   "May your Hamiltonians diagonalise in a manner truly sublime."

for data 
do
    echo "Generating Js from CGMD datafile ${data}"

awk < ${data} '{ print $1, $2, $3 }'  | ./cg-pcbm-j > ${data%.*}.edges

python DBTW.py 1000 ${data%.*}.edges

done

