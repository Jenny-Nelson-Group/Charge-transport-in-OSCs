#! /bin/bash

method="b3lyp/6-311+g*"
namemol1=$1
namemol2=$2
namecalc=$3

echo "%nproc=1
%mem=100Mb
#p punch=mo METHOD nosymm

autogen

0 1" > header


mkdir part1 
mkdir part2
mkdir dim

echo "usage xyzfile1 $namemol1 zyzfile2  $namemol2 namecala $namecalc internal: method $method"

echo "NB: at the moment you cannot specify chkpoint change by hand if needed, also _you_ gotta run the g03 jobs"

sed "s:METHOD:${method}:" header > part1/${namecalc}part1.com
cat $namemol1 >> part1/${namecalc}part1.com
awk '{printf "%s-Bq %f %f %f \n", $1, $2, $3,$4}' $namemol2 >> part1/${namecalc}part1.com
echo >>  part1/${namecalc}part1.com


sed "s:METHOD:${method}:" header > part2/${namecalc}part2.com
awk '{printf "%s-Bq %f %f %f \n", $1, $2, $3,$4}' $namemol1 >> part2/${namecalc}part2.com
cat $namemol2 >> part2/${namecalc}part2.com
echo >>  part2/${namecalc}part2.com


echo "%chk=dimer.chk" > dim/${namecalc}dim.com
sed "s:METHOD:${method}:" header >> dim/${namecalc}dim.com
cat $namemol1 >> dim/${namecalc}dim.com
cat $namemol2 >> dim/${namecalc}dim.com
echo >>  dim/${namecalc}dim.com
echo "--Link1--
%chk=dimer.chk
%mem=100Mb
#p geom(allcheck) guess(read,only) IOp(3/33=1) ${method} nosymm

" >> dim/${namecalc}dim.com



