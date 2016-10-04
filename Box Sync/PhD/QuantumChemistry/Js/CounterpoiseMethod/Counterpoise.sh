#! /bin/bash

namemol1=$1
namemol2=$2
namecalc=$3

cat HEADER_Js_CP > ${namecalc}A.com
cat $namemol1 >> ${namecalc}A.com
awk '{printf "%s-Bq %f %f %f \n", $1, $2, $3,$4}' $namemol2 >> ${namecalc}A.com
echo >>  ${namecalc}A.com


cat HEADER_Js_CP > ${namecalc}B.com
awk '{printf "%s-Bq %f %f %f \n", $1, $2, $3,$4}' $namemol1 >> ${namecalc}B.com
cat $namemol2 >> ${namecalc}B.com
echo >>  ${namecalc}B.com

echo "%chk=dimer.chk" > ${namecalc}AB.com
cat HEADER_Js_CP >> ${namecalc}AB.com
cat $namemol1 >> ${namecalc}AB.com
cat $namemol2 >> ${namecalc}AB.com
echo >> ${namecalc}AB.com
cat HEADER_Js_CP_pair >> ${namecalc}AB.com
echo >>  ${namecalc}AB.com

