for dir in C60 #M B T Tetra Penta Hexa
do

cd "${dir}"

for num in 000 #{0..9} 0{10..99} 100
do
echo $num 
for pdb in ${dir}_SimulAnneal_${num}.pdb
do
for J_cutoff in 0.001 0.005 0.01 0.05 0.1
do
 cell=` head "${pdb}" | grep CRYST1 | awk '{print $2,$3,$4}' ` # extract periodic cell vectors
 grep "^ATOM" "${pdb}" | grep " C " | awk '{print $6,$7,$8}' > "${pdb}".xyz #all the 'C'60 pseudo-atoms
 sites=` cat "${pdb}".xyz | wc -l ` #cat'ing so just have single field of lines

python ../Percolation.py ${sites} "${pdb}".xyz ${cell} ${J_cutoff} > Percolation_output_"${dir}"_"${num}"_"${J_cutoff}".txt

done
done
done

cd -

done
