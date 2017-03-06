
#!/usr/bin/env python

# Transform molcules to equivalent positions and then translate, then find molecules in the unit cell to make xyz files of each in supercell for calculation of J


import numpy as np
import sys
import Spacegroup


# --------------------------------- Read in data --------------------------------------- #

data = sys.argv[1]

# Extract raw data from xyz file (i.e. discard atom symbols)
LABELS = np.genfromtxt(data,usecols=0,dtype=str)
COORDFILE = np.genfromtxt(data)[:,1:]

# print "Coordinates: ", coordfile
N=len(COORDFILE)

# Cell lengths
a=float(sys.argv[2])
b=float(sys.argv[3])
c=float(sys.argv[4])

# Angles
alpha=np.radians(float(sys.argv[5]))
beta=np.radians(float(sys.argv[6]))
gamma=np.radians(float(sys.argv[7]))

# Volume of parallelepiped defining unit cell
omega=a*b*c*np.sqrt(1-(np.cos(alpha))**2-(np.cos(beta))**2-(np.cos(gamma))**2+2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
 
# Read in space group: Pbca, Cc, P21c, P21, Pna21, Pbcn, P212121, C2c, P_1,Pca21, C2
SG = sys.argv[8]

# File name for output files
filename = sys.argv[9]

# Conversions between fractional and real coordinates (inc for non orthogonal cells)

xyz_to_frac = [[1./a,-np.cos(gamma)/(a*np.sin(gamma)),b*c*(np.cos(alpha)*np.cos(gamma)-np.cos(beta))/(omega*np.sin(gamma))],
[0.,1./(b*np.sin(gamma)),a*c*(np.cos(beta)*np.cos(gamma)-np.cos(alpha))/(omega*np.sin(gamma))],
[0.,0.,a*b*np.sin(gamma)/omega] ]

frac_to_xyz = [  [a,b*np.cos(gamma),c*np.cos(beta)],
               [0.,b*np.sin(gamma),c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)],
               [0.,0.,omega/(a*b*np.sin(gamma))] ]


# -------------------------------------------------------------------------------------- #

def translations():

# Translation matrices

	Translations=[ [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
              [[1,0,0,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
              [[1,0,0,0],[0,1,0,1],[0,0,1,0],[0,0,0,1]],
              [[1,0,0,0],[0,1,0,0],[0,0,1,1],[0,0,0,1]],
              [[1,0,0,-1],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
              [[1,0,0,0],[0,1,0,-1],[0,0,1,0],[0,0,0,1]],
              [[1,0,0,0],[0,1,0,0],[0,0,1,-1],[0,0,0,1]],
              [[1,0,0,1],[0,1,0,1],[0,0,1,0],[0,0,0,1]],
              [[1,0,0,1],[0,1,0,-1],[0,0,1,0],[0,0,0,1]],
              [[1,0,0,-1],[0,1,0,1],[0,0,1,0],[0,0,0,1]],
              [[1,0,0,-1],[0,1,0,-1],[0,0,1,0],[0,0,0,1]],
              [[1,0,0,1],[0,1,0,0],[0,0,1,1],[0,0,0,1]],
              [[1,0,0,1],[0,1,0,0],[0,0,1,-1],[0,0,0,1]],
              [[1,0,0,-1],[0,1,0,0],[0,0,1,1],[0,0,0,1]],
              [[1,0,0,-1],[0,1,0,0],[0,0,1,-1],[0,0,0,1]],
              [[1,0,0,0],[0,1,0,1],[0,0,1,1],[0,0,0,1]],
              [[1,0,0,0],[0,1,0,1],[0,0,1,-1],[0,0,0,1]],
              [[1,0,0,0],[0,1,0,-1],[0,0,1,1],[0,0,0,1]],
              [[1,0,0,0],[0,1,0,-1],[0,0,1,-1],[0,0,0,1]],
              [[1,0,0,1],[0,1,0,-1],[0,0,1,-1],[0,0,0,1]],
              [[1,0,0,-1],[0,1,0,-1],[0,0,1,1],[0,0,0,1]],
              [[1,0,0,-1],[0,1,0,1],[0,0,1,-1],[0,0,0,1]],
              [[1,0,0,-1],[0,1,0,1],[0,0,1,1],[0,0,0,1]],
              [[1,0,0,1],[0,1,0,-1],[0,0,1,1],[0,0,0,1]],
              [[1,0,0,1],[0,1,0,1],[0,0,1,-1],[0,0,0,1]],
              [[1,0,0,1],[0,1,0,1],[0,0,1,1],[0,0,0,1]],
              [[1,0,0,-1],[0,1,0,-1],[0,0,1,-1],[0,0,0,1]]
              ]

	return Translations

# -------------------------------------------------------------------------------------- #

def transform_and_translate(coordfile):

	# Transform to fractional coordinates
	coordfile=np.inner(xyz_to_frac,coordfile).T

	#print "Fractional coordinates: ", coordfile

	#Add a column of ones for matrix multiplication as in Vesta
	coordfile=np.c_[coordfile,np.ones(N)]

	# Transform matrices

	Transforms = np.array(Spacegroup.choose(SG))

    #print Transforms,np.shape(Transforms)
    
	m=np.shape(Transforms)[0]

	all_transformations = np.zeros((1,N,3))

	for displace_n,displace in enumerate(translations()):
	    for transform_n,transform in enumerate(Transforms):
	        # Transform coordinates
			trans=np.inner(transform,coordfile).T
	        # Translate in one of 6 directions
			trans=np.inner(displace,trans).T
	        # Remove column of ones
			trans=trans[:,[0,1,2]]
	        # Scale up to real coordinates
			trans=np.inner(frac_to_xyz,trans).T
			all_transformations=np.vstack((all_transformations,trans[None,...]))
    
	return all_transformations[1:],m

# -------------------------------------------------------------------------------------- #

def check_in_cell(coordfile):
	coordfile_fractional=np.inner(xyz_to_frac,coordfile).T
	atoms=[]
    
    # Check how many atoms of each molecule are in the unit cell
    
	for i in range(0,N):
		if 0<coordfile_fractional[i,0]<1 and 0<coordfile_fractional[i,1]<1 and 0<coordfile_fractional[i,2]<1:
			atoms=np.append(atoms,i)

	return len(atoms)

# -------------------------------------------------------------------------------------- #

def unit_cell(all_molecules):
	in_cell=[0,0]
    
    # Check if any atoms in molecule are in the unit cell
    
	for i in range(0,len(all_molecules)):
		check=check_in_cell(all_molecules[i,:,:])
		if check>0:
			in_cell=np.vstack((in_cell,[i,check]))
	in_cell=in_cell[1:]

    # Sort molecules by amount in unit cell

	sorted_idx=np.argsort(in_cell[:,1])
	in_cell_sorted=in_cell[sorted_idx][::-1]

	unitcell_idx=[0,0]

    # Define one molecule in each orientation as a member of the unit cell

	for i in range(len(in_cell_sorted)):
		for j in range(1,M):
			if in_cell_sorted[i,0]%M==j and (j in unitcell_idx) == False:
				unitcell_idx=np.vstack((unitcell_idx,[j,in_cell_sorted[i,0]]))


	print unitcell_idx

	unitcell=all_molecules[unitcell_idx[:,1],:,:]

    # Outputs the molecule coordinates of all molecules in the unit cell

	return unitcell

# -------------------------------------------------------------------------------------- #

def make_supercell(UNIT_CELL):

    # Make a supercell starting with molecules from unit cell
    
    all_molecules = np.zeros((1,N,3))

    for displace_n,displace in enumerate(translations()):
    	for i in range(0,M):
            trans=np.inner(xyz_to_frac,UNIT_CELL[i,:,:]).T
            trans=np.c_[trans,np.ones(N)]
            trans=np.inner(displace,trans).T
            trans=trans[:,[0,1,2]]
            trans=np.inner(frac_to_xyz,trans).T
            all_molecules=np.vstack((all_molecules,trans[None,...]))

    return all_molecules[1:]

# -------------------------------------------------------------------------------------- #

def save_files(file_to_save,labels,FILENAME,i):
	# Put labels back to save in xyz format
	trans=np.array(zip(labels,file_to_save[:,0],file_to_save[:,1],file_to_save[:,2]),dtype=[('labels','S8'),('trans[:,0]',float),('trans[:,1]',float),('trans[:,2]',float)])
	# Save
	np.savetxt("%s_%d.xyz"%(FILENAME,i),trans,delimiter=" 	",fmt=["%s"]+["%f"]+["%f"]+["%f"])


# -------------------------------------------------------------------------------------- #


all_trans,M = transform_and_translate(COORDFILE)
unitcell=unit_cell(all_trans)
all_mols=make_supercell(unitcell)

for i in range(0,27*M):
save_files(all_mols[i,:,:],LABELS,filename,i)
