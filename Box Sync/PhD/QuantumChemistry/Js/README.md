# Calculate J

Calculating transfer integrals with Gaussian and the projective method.

Adapted from code by J.Kirkpatrick to deal with GAMESS-UK output; joejk2 12/09/07, with edits from Jarvist Frost and Beth Rice

Files required:

- MoleculeA.xyz
- MoleculeB.xyz
- HEADER_opt
- HEADER_Js
- rewrite_S_phi_E.cpp
- components.h
- get_J_interactive.m
- jkp_extract_geom.awk

Additional files for running on cx1 at Imperial College:

- launch_coms.sh
- launch_coms_special.sh

Programs required:

- Gaussian 
- Octave


1. Optimise MoleculeA and MoleculeB


cat HEADER_opt MoleculeX > MoleculeX_opt.com

echo -en "\n" >> MoleculeX_opt.com            (Gaussian requires an empty line at the end of the input file)

HEADER_opt peforms a geometry optimisation with the B3LYP level of theory and 6-31g* basis set.

Run this job with Gaussian: on cx1 this can be done with

bash launch_coms.sh MoleculeX_opt.com


2. Extract optimised geometry for each and create 3 Gaussian input files for calculating J

awk -f jkp_extract_geom.awk MoleculeX_opt.log > MoleculeX_opt.xyz

(1) cat HEADER_Js MoleculeA_opt.xyz > MoleculeA_J.com
(2) cat HEADER_Js MoleculeB_opt.xyz > MoleculeB_J.com
(3) cat HEADER_Js MoleculeA_opt.xyz MoleculeB_opt.xyz > MoleculeA_MoleculeB_J.com

Ensure the files contain a blank final line with

echo -en "\n" >> MoleculeX_J.com

N.B. HEADER_Js contains the line 

# b3lyp nosymm punch=mo iop(3/33=1)
	nosymm is necessary to prevent reorientation of molecule to standard orientation.
	punch=mo puts orbital coefficients in fort.7.  This must be moved to xxx_x.pun 
	iop(3/33=1) forces output of the overlap matrix in the xxx_x.log file
	

3. Run each of these Gaussian jobs

On cx1 this can be done with 

bash launch_coms_special.sh MoleculeX_J.com


4. Compile and run 'rewrite_S_Phi_E.x' with the prefix to gaussian calculation (1), (2)
and (3) !! IN THIS ORDER !!.  This program expects to find 3 log files and 3 pun files.

g++ rewrite_S_Phi_E.cpp -o rewrite_S_Phi_E.x (requires the header file components.h)


5. You should now have 7 output .txt files:
	Evls_pair.txt  MOs_1.txt  MOs_2.txt  MOs_pair.txt  S_1.txt  S_2.txt
	S_pair.txt
	

6. Run 'octave get_J_interactive.m'.  You will be prompted for the number of MOs and
the HOMO of a monomer, which can be found in the log file under 'basis functions,' and 'alpha electrons' respectively.


NOTE:

The script file Calculate_J.sh automates steps 4, 5 and 6 (given rewrite_S_Phi_E is already compiled), while the script file Extract_Js.sh can be used to extract the J for HOMO or LUMO coupling as desired. 





