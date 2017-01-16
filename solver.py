#####################################################################
#
# solver.py 
# This script solves for the correct conformal frame for the
# handlebody geometries. The input is the boundary circles and the
# output is the solution. 
#
#####################################################################

import subprocess
import numpy as np

# Import variables from makemesh function
mesh_out = subprocess.check_output(['./makemesh','0','0','-1']).split("break")

# Import the coordinates
coors = np.fromstring(mesh_out[0][1:-2], sep=", ")
coors = coors.reshape((coors.size/2,2))

print coors[-1]

