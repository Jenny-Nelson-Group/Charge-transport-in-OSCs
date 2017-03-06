# Calculate Mobility

The python code Mobility.py calculates the charge carrier mobility in a crystalline system where the Marcus equation for rates applies (J<<lambda). The code finds the angular dependent mobility in the ab, bc and ac planes, and plots on a polar plot.

The following inputs are required:

1. Locationsfile of xyz coordinate of centres of molecules in crystal

2. Jfile of electronic couplings between molecules, given in rows of $molecule1 $molecule2 $J

3. Magnitude of electric field applied, in A/cm

4. Filename for saving plots