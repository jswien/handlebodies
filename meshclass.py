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
        self.mesh_length = self.coors.size/2
        self.elements_to_nodes = self._build_elements_to_nodes(makemesh_output)
        self.nodes_to_elements = self._build_nodes_to_elements(makemesh_output)
        self.bndy_elements = self._build_bndy_elements(makemesh_output)
        self.neighborhoods = self._build_neighborhoods()
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

    def _build_elements_to_nodes(self, makemesh_output):
        """Build list of nodes contained in each element"""
        elements_to_nodes = np.fromstring(makemesh_output[1][2:-2], dtype=int,sep=",")-1
        return  elements_to_nodes.reshape(elements_to_nodes.size/6,6).tolist()

    def _build_nodes_to_elements(self, makemesh_output):
        """Build the array of elements connected to each node"""
        connection_list =  makemesh_output[3][3:-3].split("}, {")
        nodes_to_elements = []
        for connection in connection_list:
            nodes_to_elements.append([int(n)-1 for n in connection.split(",")])
        return nodes_to_elements

    def _build_bndy_elements(self, makemesh_output):
        """Read in list of boundary elements"""
        return  np.fromstring(makemesh_output[2][2:-2], dtype=int,sep=",")

    def _build_neighborhoods(self):
        """Build a list of nodes adjacent to node i"""
        neighborhoods = [[]]*self.mesh_length
        for node in xrange(self.mesh_length):
            neighborhoods[node]= [self.elements_to_nodes[elem] for elem in self.nodes_to_elements[node]]
            neighborhoods[node]= set([elem for block in neighborhoods[node] for elem in block])
        return neighborhoods
            

    def _build_psi_metric(self):
        """Build matrix of integrals psi_i psi_j"""
        pass


mesh = Mesh(np.array([[1,0,0],[0,1,0],[0,0,-1]]), np.array([0.005,4]))

print mesh.neighborhoods[9]
