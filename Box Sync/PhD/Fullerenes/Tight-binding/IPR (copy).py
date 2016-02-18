#!/usr/bin/env python


import numpy as np
import sys

F_mag=1.0

for theta in xrange(0,360,10):
	f1[theta]=F_mag*np.cos(theta)





# Save fields
np.savetxt("Field_ab.txt",Field,delimiter=' ',newline='\n')








