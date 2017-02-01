#! /bin/bash

file=$1
namecalc=$2

awk -v len=$(wc -l < "${file}") -v fileA="${file%.*}"A -v fileB="${file%.*}"B '
  BEGIN    { outfile = fileA  } 
  NR>len/2 { outfile = fileB  }
           { print $0 > outfile }
' "${file}"

echo "File split"

cat ../HEADER_Js_CP_bonded > ${namecalc}A.com
cat "${file%.*}"A  >> ${namecalc}A.com
awk '{printf "%s-Bq %f %f %f \n", $1, $2, $3,$4}' "${file%.*}"B >> ${namecalc}A.com
echo >>  ${namecalc}A.com

echo "First com file made"

cat ../HEADER_Js_CP_bonded > ${namecalc}B.com
awk '{printf "%s-Bq %f %f %f \n", $1, $2, $3,$4}' "${file%.*}"A  >> ${namecalc}B.com
cat "${file%.*}"B  >> ${namecalc}B.com
echo >>  ${namecalc}B.com

echo "Second com file made"

echo "%chk=dimer.chk" > ${namecalc}AB.com
cat ../HEADER_Js_CP >> ${namecalc}AB.com
cat "${file%.*}"A  >> ${namecalc}AB.com
cat "${file%.*}"B  >> ${namecalc}AB.com
echo >> ${namecalc}AB.com
cat ../HEADER_Js_CP_pair >> ${namecalc}AB.com
echo >>  ${namecalc}AB.com

echo "Done"
