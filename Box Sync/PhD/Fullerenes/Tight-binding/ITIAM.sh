for pdb
do
 cell=` head "${pdb}" | grep CRYST1 | awk '{print $2,$3,$4}' ` # extract periodic cell vectors
 grep "^ATOM" "${pdb}" | grep " C " | awk '{print $6,$7,$8}' > "${pdb}".xyz #all the 'C'60 pseudo-atoms
 sites=` cat "${pdb}".xyz | wc -l ` #cat'ing so just have single field of lines


for alpha in 0.32
do
for dx in 0
do


python ITIAM.py ${sites} "${pdb}".xyz ${cell} $alpha $dx

done
done



#for alpha in 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 02
#do
#for dx in 0 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 0.011 0.012 0.013 0.014 0.015 0.016 0.017 0.018 0.019 0.02 0.021 0.022 0.023 0.024 0.025
#do
#python ITIAM.py ${sites} tmp.xyz ${cell} $alpha $dx > outputs_${alpha}_${dx}.dat
#grep "Localised" outputs_${alpha}_${dx}.dat > localise_${alpha}_${dx}.dat
#rm outputs_${alpha}_${dx}.dat
#if [ -s localise_${alpha}_${dx}.dat ]
# then
#break
# else
#rm localise_${alpha}_${dx}.dat
#fi

#done
#done


done



