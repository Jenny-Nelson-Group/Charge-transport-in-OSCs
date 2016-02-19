for dir in Penta #C60 M B T Tetra Penta Hexa
do

cd "${dir}"

for num in 0{75..99} 100 #00{0..9} 0{10..99} 100
do

for pdb in ${dir}_SimulAnneal_${num}.pdb
do

echo $num

 cell=` head "${pdb}" | grep CRYST1 | awk '{print $2,$3,$4}' ` # extract periodic cell vectors
 grep "^ATOM" "${pdb}" | grep " C " | awk '{print $6,$7,$8}' > "${pdb}".xyz #all the 'C'60 pseudo-atoms
 sites=` cat "${pdb}".xyz | wc -l ` #cat'ing so just have single field of lines

evals=$(python ../Evals.py ${sites} "${pdb}".xyz ${cell} 0)

echo $evals>>Eigenvalues_E0.dat


done
done

cd -

done



