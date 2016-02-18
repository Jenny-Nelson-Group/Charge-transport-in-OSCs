for dx in 0.0
do
for alpha in 0.2
do
for scfstep in {1..10}
do

python ITIAM.py 10 10_line.xyz 10 1 1 $alpha $dx $scfstep

done
done
done





