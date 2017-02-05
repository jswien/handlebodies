#####################################################################
#
# phi_solver.py
#
# Solves for the constant curvature metric for the Riemann surface
# defined by a set of circles with identifications.
#
#####################################################################

import time

import numpy as np
from numpy import pi as pi
from numpy import sqrt as sqrt

from mesh_class import Mesh as Mesh

print time.ctime()
mesh = Mesh(np.array([[1,0,0],[0,1,0],[0,0,-1]]), np.array([0.005,4]))
print time.ctime()
print len(mesh.coors)

print np.ones(mesh.mesh_length).dot(mesh.psi_metric.dot(np.ones(mesh.mesh_length))) - pi/4

