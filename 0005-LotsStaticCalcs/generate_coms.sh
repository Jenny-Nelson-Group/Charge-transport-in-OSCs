# generate some Efields for Gaussian calcs.

for i
do

for FIELD in \
    "z+0" "z+1" "z+3" "z+10" "z+32" "z+100" "z+316" "z+1000" "z+3162" "z+9999" \
    "z-1" "z-3" "z-10" "z-32" "z-100" "z-316" "z-1000" "z-3162" "z-9999" \
    "y+0" "y+1" "y+3" "y+10" "y+32" "y+100" "y+316" "y+1000" "y+3162" "y+9999" \
    "y-1" "y-3" "y-10" "y-32" "y-100" "y-316" "y-1000" "y-3162" "y-9999" \
    "x+0" "x+1" "x+3" "x+10" "x+32" "x+100" "x+316" "x+1000" "x+3162" "x+9999" \
    "x-1" "x-3" "x-10" "x-32" "x-100" "x-316" "x-1000" "x-3162" "x-9999" 


    #"x+0" "x+1" "x+10" "x+100" "x+1000" "y+1" "y+10" "y+100" "y+1000"
do
 filename=` echo "${i%.*}"_${FIELD}.com | sed -e s/+/p/ -e s/-/m/ `
 echo $filename

 cat HEADER | sed "s/FIELD/${FIELD}/" > $filename
 cat "${i}" >> $filename
 echo "" >> $filename
done
done
