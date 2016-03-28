# Calculate Total Energy, HOMOs and LUMOs

These can be calculated with Gaussian by running an optimisation (see HEADER_opt for commands) and extracting the total energy and Kohn-Sham orbitals. However, the LUMO from this is not correct and a better method is the delta SCF method (see this by Jarvist Frost for more info http://jarvist.github.io/post/2016-03-17-delta-scf/ )

Gaussian single point energy calculations are done on the neutral geometry molecule for the neutrally charged molecule, cation and anion. Use HEADER_neutral, HEADER_cation and HEADER_anion. 

Then the script Extract_values.sh can be used where Type = HOMO or LUMO means the Kohn-Sham HOMO and LUMOs, and HOMO_SCF or LUMO_SCF means the HOMO and LUMO with the delta SCF method (all in eV). Type= Energy extracts the total energy in Hartrees.





