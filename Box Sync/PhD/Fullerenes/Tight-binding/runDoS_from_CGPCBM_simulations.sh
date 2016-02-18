#!/bin/sh
#  runDoS_from_CGPCBM_simulations.sh - generates Js (.edges) file with ./cg-pcbm-j for specified PDB PCBM CG files
#       & runs python TB ToS calculator (DBTW.py)
#   "May your Hamiltonians diagonalise in a manner truly sublime."

for data 
do
    echo "Generating Js from CGMD datafile ${data}"
#OK, let's generate our edges (connection data)
# # with our dirty little C code program, and a bit of AWK manipulation
pdbcat -fields "${data}"  | grep " C " | awk 'BEGIN{print 1000, 30}{print $10,$11,$12}'  | ./cg-pcbm-j > ${data%.*}.edges

python DBTW.py 1000 ${data%.*}.edges

done

