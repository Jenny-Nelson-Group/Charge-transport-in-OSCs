# generate some Efields for Gaussian calcs.

for i
do

for FIELD in "x+0" "x+1" "x+10" "x+100" "x+1000" "y+1" "y+10" "y+100" "y+1000"
do
 filename=` echo "${i%.*}"_${FIELD}.com | sed -e s/+/p/ -e s/-/m/ `
 echo $filename

 cat HEADER | sed "s/FIELD/${FIELD}/" > $filename
 cat "${i}" >> $filename
 echo "" >> $filename
done
done
