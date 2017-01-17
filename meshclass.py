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

        self.psi_metric=self._build_global_metric(local_psi_metric)
        self.Npsi_metric=self._build_global_metric(local_Npsi_metric)

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
        elements_to_nodes = np.fromstring(makemesh_output[1][2:-2], 
                                          dtype=int,sep=",")-1
        return  elements_to_nodes.reshape(elements_to_nodes.size/6,6).tolist()

    def _build_nodes_to_elements(self, makemesh_output):
        """Build the array of elements connected to each node"""
        connection_list =  makemesh_output[3][3:-3].split("}, {")
        nodes_to_elements = []
        for connection in connection_list:
            nodes_to_elements.append([int(n)-1 for n 
                                      in connection.split(",")])
        return nodes_to_elements

    def _build_bndy_elements(self, makemesh_output):
        """Read in list of boundary elements"""
        return  np.fromstring(makemesh_output[2][2:-2],
                              dtype=int,sep=",")

    def _build_neighborhoods(self):
        """Build a list of nodes adjacent to each node"""
        neighborhoods = [[]]*self.mesh_length
        for node in xrange(self.mesh_length):
            neighborhoods[node]= [
                self.elements_to_nodes[elem] 
                for elem in self.nodes_to_elements[node]
                ]
            neighborhoods[node]= set([elem for 
                                      block in neighborhoods[node]for
                                      elem in block])
        return neighborhoods
    
    
    def _build_global_metric(self, metric_function):
        """Build global metric from local function"""
        metric = np.zeros((self.mesh_length,self.mesh_length))
        for n_i in xrange(self.mesh_length):
            for n_j in self.neighborhoods[n_i]:
                for elem in list(set(self.nodes_to_elements[n_i])&
                                 set(self.nodes_to_elements[n_j])
                                 ):
                    n_i_pos = self.elements_to_nodes[elem].index(n_i)
                    n_j_pos = self.elements_to_nodes[elem].index(n_j)
                    metric[n_i,n_j]+=metric_function(
                        np.array([self.coors[node] for node 
                                  in self.elements_to_nodes[elem][:3]
                                 ])
                        )[n_i_pos,n_j_pos]
        return metric
                               
def local_psi_metric(vertices):
    jac = np.absolute(vertices[0,0]*(vertices[1,1]-vertices[2,1])
                      +vertices[1,0]*(vertices[2,1]-vertices[0,1])
                      +vertices[2,0]*(vertices[0,1]-vertices[1,1])
                      )
    vv = jac/60
    vvp = -jac/360
    vn = -jac/90
    nn = 4*jac/45
    nnp = 2*jac/45
    return np.array([[vv,vvp,vvp,0,vn,0],
                     [vvp,vv,vvp,0,0,vn],
                     [vvp,vvp,vv,vn,0,0],
                     [0,0,vn,nn,nnp,nnp],
                     [vn,0,0,nnp,nn,nnp],
                     [0,vn,0,nnp,nnp,nn],
                     ])

def local_Npsi_metric(vertices):
    jac = np.absolute(vertices[0,0]*(vertices[1,1]-vertices[2,1])
                      +vertices[1,0]*(vertices[2,1]-vertices[0,1])
                      +vertices[2,0]*(vertices[0,1]-vertices[1,1])
                      )
    area = 4*((vertices[0,0]-vertices[1,0])**2
              +(vertices[0,0]-vertices[2,0])**2
              +(vertices[2,0]-vertices[1,0])**2
              +(vertices[0,1]-vertices[1,1])**2
              +(vertices[0,1]-vertices[2,1])**2
              +(vertices[2,1]-vertices[1,1])**2
             )
    c00 = ((vertices[1,0]-vertices[2,0])**2
           +(vertices[1,1]-vertices[2,1])**2)
    c11 = ((vertices[0,0]-vertices[2,0])**2
           +(vertices[0,1]-vertices[2,1])**2)
    c22 = ((vertices[0,0]-vertices[1,0])**2
           +(vertices[0,1]-vertices[1,1])**2)
    c01 = ((vertices[0,0]-vertices[2,0])
           *(vertices[1,0]-vertices[2,0])
           +(vertices[0,1]-vertices[2,1])
           *(vertices[1,1]-vertices[2,1]))
    c02 = ((vertices[0,0]-vertices[1,0])
           *(vertices[2,0]-vertices[1,0])
           +(vertices[0,1]-vertices[1,1])
           *(vertices[2,1]-vertices[1,1]))
    c12 = ((vertices[1,0]-vertices[0,0])
           *(vertices[2,0]-vertices[0,0])
           +(vertices[1,1]-vertices[0,1])
           *(vertices[2,1]-vertices[0,1]))
    
    return 1/jac*1/6*np.array([[3*c00,c01,c02,-4*c01,0,-4*c02],
                               [c01,3*c11,c12,-4*c01,-4*c12,0],
                               [c02,c12,3*c22,0,-4*c12,-4*c02],
                               [-4*c01,-4*c01,0,area,-8*c02,-8*c12],
                               [0,-4*c12,-4*c12,-8*c02,area,-8*c01],
                               [-4*c02,0,-4*c02,-8*c12,-8*c01,area],
                               ])

mesh = Mesh(np.array([[1,0,0],[0,1,0],[0,0,-1]]), np.array([0.005,4]))

print np.ones(mesh.mesh_length).dot(mesh.psi_metric.dot(np.ones(mesh.mesh_length)))-np.pi/4

print np.ones(mesh.mesh_length).dot(mesh.Npsi_metric.dot(np.ones(mesh.mesh_length)))
