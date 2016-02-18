#!/bin/sh

Name=$1
Locationsfile=$2
Jfile=$3
Fieldfile=$4
cell_a=$5
cell_b=$6
cell_c=$7
Field_mag=$8
FIELD=$9

Angle=Angles

for field_line in {1..36}
do

field=$(sed -n "$field_line"p $Fieldfile) 
angle=$(sed -n "$field_line"p $Angle)

echo $field

Mobility=$(python Mobility.py 108 $Jfile $Locationsfile $cell_a $cell_b $cell_c $Field_mag $field $angle) 

echo $Mobility

echo $Mobility >> "${Name}"_mobility_vs_angle_"${FIELD}".dat
done
 
#rm locations.txt



#for i in {0..107}
#do
#location=$(python Find_centres.py ./R2/R2_"$i".xyz 41)
#echo "$location">>locations_R2.txt
#done

#Fieldfile=Field_ac_rads.txt
#Jfile=./Js/Js_HOMO_R2_cell.txt
#Angle=Angles

#for field_line in {1..36}
#do

#field=$(sed -n "$field_line"p $Fieldfile) 
#angle=$(sed -n "$field_line"p $Angle)

#echo $field

#Mobility=$(python Mobility.py 107 $Jfile locations_R2.txt 9.1753 10.3659 17.7857 $field $angle) 

#echo $Mobility

#echo $Mobility >> R2_mobility_vs_angle_acplane.txt
#echo -en "\n" >> EP13_mobility_vs_angle_abplane.txt

#done
 
#rm locations.txt


