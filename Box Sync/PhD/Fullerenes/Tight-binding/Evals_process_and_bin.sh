for dir in C60 M B T Tetra Penta Hexa
do

cd "${dir}"

sed 's/ /\n/g' Eigenvalues_E0.dat > Evals_"${dir}".dat

python ../Bin_data.py Evals_"${dir}".dat "${dir}" 

cd -
done





