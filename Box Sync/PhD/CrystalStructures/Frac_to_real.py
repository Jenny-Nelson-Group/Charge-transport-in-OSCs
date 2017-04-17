#!/usr/bin/env python

# Convert fractional coordintaes to real (inc. for non orthogonal unit cells)

import numpy as np
import sys

data = sys.argv[1]

# Extract data from xyz file 
coordfile = np.loadtxt(data)

# print "Coordinates: ", coordfile=
n=len(coordfile)

# Cell lengths
a=float(sys.argv[2])
b=float(sys.argv[3])
c=float(sys.argv[4])

# Angles
alpha=np.radians(float(sys.argv[5]))
beta=np.radians(float(sys.argv[6]))
gamma=np.radians(float(sys.argv[7]))

# Volume of parallelepiped defining unit cell
omega=a*b*c*np.sqrt(1-np.cos(alpha)**2-np.cos(beta)**2-np.cos(gamma)**2+2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))

# File name for output files
filename = sys.argv[8]

# Matrix to convert to xyz coordinates from fractional (inc. for non orthogonal)
frac_to_xyz = [  [a,b*np.cos(gamma),c*np.cos(beta),],
               [0,b*np.sin(gamma),(c*np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)],
               [0,0,omega/(a*b*np.sin(gamma))] ]

# Transform to real coordinates
coordfile=np.inner(frac_to_xyz,coordfile).T

# Define atom symbols for hexahelicene
labels=['C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C',
'C','C','C','C','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H']


# Make xyz file with atom symbols
xyz_file=np.array(zip(labels,coordfile[:,0],coordfile[:,1],coordfile[:,2]),dtype=[('labels','S8'),('coordfile[:,0]',float),('coordfile[:,1]',float),('coordfile[:,2]',float)])

# Save
np.savetxt("%s"%(filename),xyz_file,delimiter=" ",fmt=["%s"]+["%f"]+["%f"]+["%f"])


