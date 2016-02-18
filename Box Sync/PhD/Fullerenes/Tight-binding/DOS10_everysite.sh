for dx in 0.01
do
for alpha in 0.25
do
for state in 0
do

python ITIAM_everysite.py 10 10_line.xyz 10 1 1 $alpha $dx $state

done
done
done





