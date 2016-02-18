#!/bin/sh
#  runDoS_from_CGPCBM_simulations.sh - generates Js (.edges) file with ./cg-pcbm-j for specified PDB PCBM CG files
#       & runs python TB ToS calculator (DBTW.py)
#   "May your Hamiltonians diagonalise in a manner truly sublime."

for data 
do
cell=` head "${data}" | grep CRYST1 | awk '{print $2}' ` # extract periodic cell length
#echo "Generating Js from CGMD datafile ${data}"
#OK, let's generate our edges (connection data)
# # with our dirty little C code program, and a bit of AWK manipulation
#awk -v cell="$cell"


./pdbcat -fields "${data}"  | grep " C " |  awk 'BEGIN{print '$cell', 1000, 20}{print $10,$11,$12}' | ./cg-pcbm-j > ${data%.*}.edges

./pdbcat -fields "${data}"  | grep " C " | awk '{print $10,$11,$12}'  | ./cg-pcbm-pos > ${data%.*}.pos

python DBTW.py 1000 ${data%.*}.edges ${data%.*}.pos

done