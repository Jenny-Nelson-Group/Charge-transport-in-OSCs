for dir in C60 M B
do

for J_cutoff in 0.0001 0.0005 0.001 0.005 0.01 0.02 0.03 0.04 0.05
do
for num in 00{1..9} 0{10..99} 100
do

grep "max cluster size for J_cutoff is " Percolation_output_"${dir}"_"${num}"_"${J_cutoff}".txt > temp.txt

Perc=$(grep -Eo '\<[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\>' temp.txt)

echo $Perc >> Percolation_for_"${J_cutoff}"_"${dir}".txt

done

Percolation=$(awk '{ total += $1; count++ } END { print total/count }' Percolation_for_"${J_cutoff}"_"${dir}".txt)

echo $J_cutoff $Percolation >> Percolation_vs_J_cutoff_"${dir}".txt

done

done
