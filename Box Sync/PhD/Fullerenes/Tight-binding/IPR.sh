for pdb in C60_SimulAnneal_051.pdb M_SimulAnneal_050.pdb B_SimulAnneal_050.pdb T_SimulAnneal_050.pdb tetra_box.pdb pentakis_box.pdb hexakis_box.pdb
do
 cell=` head "${pdb}" | grep CRYST1 | awk '{print $2,$3,$4}' ` # extract periodic cell vectors
 grep "^ATOM" "${pdb}" | grep " C " | awk '{print $6,$7,$8}' > "${pdb}".xyz #all the 'C'60 pseudo-atoms
 sites=` cat "${pdb}".xyz | wc -l ` #cat'ing so just have single field of lines


for LUMO in -3.7
do

python IPR.py ${sites} "${pdb}".xyz ${cell} $LUMO

done
done

