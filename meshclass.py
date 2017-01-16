#####################################################################
#
# meshclass.py
#
# Defines the class 'mesh' to be used with FEM methods.
# Input is a np array of boundary circles.
#
#####################################################################

import subprocess
import numpy as np

class Mesh(object):
    
    def __init__(self, circles, accuracy=np.array([0.005,4])):
        self.circles = circles
        self.accuracy = accuracy
        
        makemesh_output = self._call_makemesh(
            accuracy,
            self.circles.reshape(self.circles.size))
        
        self.coors = self._build_coors(makemesh_output)
        self.mesh_elements = self._build_mesh_elements(makemesh_output)
        self.bndy_elements = self._build_bndy_elements(makemesh_output)
        self.node_connectivity = self._build_connectivity(makemesh_output)
    
    def _call_makemesh(self, *arguments):
        """Run the makemesh script"""
        command_list = np.array(['./makemesh'])
        for arg in arguments:
            command_list = np.append(command_list, arg)
        return subprocess.check_output(command_list).split("break")

    def _build_coors(self, makemesh_output):
        """Read in mesh coordinates from makemesh output"""
        coors = np.fromstring(makemesh_output[0][1:-2], sep=",")
        return  coors.reshape((coors.size/2,2))

    def _build_mesh_elements(self, makemesh_output):
        """Read in list of mesh elements"""
        mesh_elements = np.fromstring(makemesh_output[1][2:-2], dtype=int,sep=",")
        return  mesh_elements.reshape(mesh_elements.size/6,6)

    def _build_bndy_elements(self, makemesh_output):
        """Read in list of boundary elements"""
        return  np.fromstring(makemesh_output[2][2:-2], dtype=int,sep=",")
    def _build_connectivity(self, makemesh_output):
        """Read in the connectivity array"""
        connection_list =  makemesh_output[3][3:-3].split("}, {")
        node_connectivity = []
        for connection in connection_list:
            node_connectivity = node_connectivity + [[int(n) for n in connection.split(",")]]
        return node_connectivity


mesh = Mesh(np.array([[1,0,0],[0,1,0],[0,0,-1]]), np.array([0.005,4]))

