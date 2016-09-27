#!/bin/sh

Name=$1
Locationsfile=$2
Jfile=$3
N=$4
Field_mag=$5

#Angle=Angles

#for field_line in {1..36}
#do

#field=$(sed -n "$field_line"p $Fieldfile) 
#angle=$(sed -n "$field_line"p $Angle)

for theta in 0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360
do
for phi in 0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360
do

#echo $theta

#echo $phi

Mobility=$(python Mobility.py $N $Locationsfile $Jfile $Field_mag $theta $phi)

echo $Mobility

echo $Mobility >> "${Name}"_mobility_vs_angle.dat
done
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


