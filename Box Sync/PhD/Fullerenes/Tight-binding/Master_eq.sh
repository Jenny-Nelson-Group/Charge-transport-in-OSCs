for pdb
do
 cell=` head "${pdb}" | grep CRYST1 | awk '{print $2,$3,$4}' ` # extract periodic cell vectors
 grep "^ATOM" "${pdb}" | grep " C " | awk '{print $6,$7,$8}' > "${pdb}".xyz #all the 'C'60 pseudo-atoms
 sites=` cat "${pdb}".xyz | wc -l ` #cat'ing so just have single field of lines

python J_"${pdb}" ${sites} LUMOs_nodis.dat ${cell} Evecs_localising_all_states_"${pdb}".xyz.dat


done


