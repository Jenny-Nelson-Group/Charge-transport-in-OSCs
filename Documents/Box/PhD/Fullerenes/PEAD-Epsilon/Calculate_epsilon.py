#!/usr/bin/env python

# Calculate the high frequency (optical) and static dielectric constants given Dipole moments from a Gaussian log file

import numpy as np
import re
import sys

all_data_field=np.zeros((3,5))
all_data_0=np.zeros((3,5))

j=0

# Read in and process log file to extract floats (and ignore 60 from C60)
with open(sys.argv[1],'r') as f:
	for line in f:
		checkline=re.search('Dipole moment',line)
		if not checkline:
			strings = re.findall(r"[-+]?\d*\.\d+|\d+",line)
			data = [float(i) for i in strings]
			print data
			if data[0]==60.0:
				data.remove(60.0)
			all_data_field[j,:]=data
			j+=1

print all_data_field

k=0

with open(sys.argv[2],'r') as f:
	for line in f:
		checkline=re.search('Dipole moment',line)
		if not checkline:
			strings = re.findall(r"[-+]?\d*\.\d+|\d+",line)
			data = [float(i) for i in strings]
			if data[0]==60.0:
				data.remove(60.0)
			all_data_0[k,:]=data
			k+=1

print all_data_0

# Define constants and conversion factors

eps0 = 8.854187817620*10**-12    # C/Vm
D = 3.33564*10**-30     # Debye in C.m

# Volume of molecule

Vol = 670.423*10**-30   # m^3

# Calculate field in V/m

Field = all_data_field[:,0]*0.0001*5.14220652*10**11        # V/m

print Field

# Make polarisation tensor

DM_field=np.diag(all_data_field[:,1:4])*D
DM_0=np.diag(all_data_0[:,1:4])*D

P_field=DM_field/Vol
P_0=DM_0/Vol

print P_field
print P_0


# Calculate epsilon

eps=np.zeros((3,3))

for i in range(0,3):
    for j in range(0,3):
		eps[i,j] = (4*np.pi/eps0)*(P_field[i]-P_0[i])/Field[j]

print "High freq eps: ", eps

print "Trace/3: ", np.trace(eps)/3
