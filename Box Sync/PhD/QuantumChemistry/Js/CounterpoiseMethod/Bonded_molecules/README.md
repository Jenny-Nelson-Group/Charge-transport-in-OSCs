# Calculate the transfer integral J between two bonded monomers 

For a simple calculation, run 

bash Bonded_counterpoise.sh $xyz_dimer $filename

where $xyz_dimer is an xyz file containing the atoms and coordinates of the dimer and $filename is what you want to call the com files 

This yields 3 Gaussian com files (A, B and AB), which can be run with launch_coms.sh or similar.

Once they have finished running, do

python Counterpoise_method.py $filename_A.log $filename_B.log $filename_AB.log



If you want to change the dihedral angle between the monomers, start with a com file of the dimer with the header

#p opt=ModRedundant b3lyp/6-31g(d) 

and with the last two lines reading  

D	$1	$2	$3	$4	=DIH.0	B
D	$1	$2	$3	$4	F

where $1, $2, $3 and $4 refer to the numbers of the four atoms of the dihedral of interest. These can be found in Gaussview.  


Now run

bash Change_dihedral_angle.sh $Dihedral_com_file

on this com file to yield com files for dihedral angels at 10 degree intervals. Run these files and then extract the xyz coordinates from the resultant log files with

bash Extract_xyz_from_log.sh $file_root


Now run 

bash Make_J_files_bonded.sh $file_root $filename

which calls Bonded_counterpoise.sh to make A, B and AB com files to run for each dihedral angle, and run these com files.









