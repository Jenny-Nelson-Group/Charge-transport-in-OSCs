#!/usr/bin/env python

# Calculate the high frequency (optical) and static dielectric constants given Dipole moments from a Gaussian log file

import numpy as np
import re
import sys

# Read in and process log file to extract floats (and ignore 60 from C60)
with open(sys.argv[1],'r') as f:
    text = f.read()
    strings = re.findall(r"[-+]?\d*\.\d+|\d+",text)
    nums = [float(i) for i in strings]

dm_field=np.array(nums[-5:])
dm_field=np.delete(dm_field,-1)

print dm_field

with open(sys.argv[2],'r') as f:
    text = f.read()
    strings = re.findall(r"[-+]?\d*\.\d+|\d+",text)
    nums = [float(i) for i in strings]

dm_0=np.array(nums[-5:])
dm_0=np.delete(dm_0,-1)

print dm_0

# Define constants and conversion factors

eps0 = 8.854187817620*10**-12    # C/Vm
D = 3.33564*10**-30     # Debye in C.m

# Volume of molecule

Vol = 670.423*10**-30   # m^3

# Calculate field in V/m
Field=np.zeros(3)

Field[0] = dm_field[0]*0.0001*5.14220652*10**11        # V/m
print Field

# Make polarisation tensor

DM_field=dm_field[1:4]*D
DM_0=dm_0[1:4]*D

P_field=DM_field/Vol
P_0=[0,0,0]

print P_field
print P_0

# Calculate epsilon

eps=np.zeros((3,3))

for i in range(0,3):
    for j in range(0,3):
        eps[i,j] = (4*np.pi/eps0)*(P_field[i]-P_0[i])/Field[j]

print eps













