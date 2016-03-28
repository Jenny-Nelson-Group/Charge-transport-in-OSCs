# Calculate Mobility

The script Get_mobility.sh calculates the charge carrier mobility in a crystalline system where the Marcus equation for rates applies (J<<lambda) with the python script Mobility.py

The following inputs are required:

1. Name for output file

2. Locationsfile of xyz coordinates of atoms in molecule

3. Jfile of electronic couplings between molecules, given in rows of $molecule1 $molecule2 $J

4. FieldFile of field strengths in 36 equally spaced directions (Given in Field_ab_rads.txt Field_bc_rads.txt Field_ac_rad.txt for the ab, bc and ac planes)

5, 6, 7. cell_a, cell_b, cell_c dimensions of unit cell, in Angstroms

8. Magnitude of electric field applied, in A/cm

9. Field name for saving