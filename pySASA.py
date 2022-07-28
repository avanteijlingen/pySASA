# -*- coding: utf-8 -*-

""" Van der Waals radii in [A] taken from:
A cartography of the van der Waals territories
S. Alvarez, Dalton Trans., 2013, 42, 8617-8636
DOI: 10.1039/C3DT50599E
"""

import MDAnalysis as mda
from ase import Atoms
import numpy as np
import pandas
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from sklearn.metrics import euclidean_distances


#########################
# VALIDATION
#Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe
#http://curie.utmb.edu/getarea.html
#PROBE 1.4 A
# Total  area/energy         =         1366.31
# Number of surface atoms    =           66
# Number of buried atoms     =           12
#########################


class sssSASA:

    def generate_sphere_points(self, n):
       """
       Returns list of coordinates on a sphere using the Golden-
       Section Spiral algorithm.
       """
       points = []
       inc = np.pi * (3 - np.sqrt(5))
       offset = 2 / float(n)
       for k in range(int(n)):
          y = k * offset - 1 + (offset / 2)
          r = np.sqrt(1 - y*y)
          phi = k * inc
          point = np.array((np.cos(phi)*r, y, np.sin(phi)*r), dtype=np.float64, copy=True)
          points.append(point)
       return points
    
    def find_neighbor_indices(self, atoms, coords, probe, k):
        """
        Returns list of indices of atoms within probe distance to atom k. 
        """
        radius = self.vdw_radii.at[atoms[k], "vdw_radius"]
        neighbor_indices = []
        d = euclidean_distances(coords[k].reshape(1,3), coords)
        for i in range(d.shape[1]):
            if i == k:
                continue
            radius_i = self.vdw_radii.at[atoms[i], "vdw_radius"]
            if d[0][i] < radius + radius_i + probe: #+probe twice?
                neighbor_indices.append(i)
        return neighbor_indices

    def VisSphere(self):
        x = np.vstack(self.sphere_points)[:,0]
        y = np.vstack(self.sphere_points)[:,1]
        z = np.vstack(self.sphere_points)[:,2]
        ax = plt.figure().add_subplot(projection='3d')
        #ax.scatter3D(x, y, z, c=z, cmap='Greens')
        ax.scatter3D(*self.sphere_points.T)
        plt.show()
        
    def calcSASA(self):
        for i in range(0, self.pos.shape[0]):
            neighbor_indices = self.find_neighbor_indices(self.atoms, self.pos, self.radius_probe, i)
            n_neighbor = len(neighbor_indices)
            j_closest_neighbor = 0
            radius = self.radius_probe + self.vdw_radii.at[self.atoms[i], "vdw_radius"]
            
            n_accessible_point = 0
            for point in self.sphere_points:
                is_accessible = True
                # Move sphere point to atomic coord and scale by atomic radius
                test_point = (point*radius) + self.pos[i]
        
                # speed up over np.arange(X, Y) as it starts at a more likely indice meaning less for loop iterations
                # i.e. instead of [0,1,2,3,4,5] it might be [3,4,5,0,1]
                cycled_indices = np.hstack((np.arange(j_closest_neighbor, n_neighbor),
                                            np.arange(0, j_closest_neighbor)))
        
                for j in cycled_indices:
                    pos_j = self.pos[neighbor_indices[j]]
                    radius_j = self.radius_probe + self.vdw_radii.at[self.atoms[j], "vdw_radius"]
                    diff = np.linalg.norm(pos_j - test_point)
                    if diff**2 < radius_j**2:
                        j_closest_neighbor = j
                        is_accessible = False
                        break
                if is_accessible:
                    n_accessible_point += 1
                    accessible_points = np.vstack((self.accessible_points, test_point))
                
            area = self.area_per_point*n_accessible_point*radius**2 
            self.areas.loc[i] = [area, self.atoms[i], self.mol[i].segid, 
                                 self.mol[i].resname, self.mol[i].resid, 
                                 self.vdw_radii.at[self.atoms[i], "vdw_radius"]]
    
    def writeConnolly(self, fname="ConnollySurface.xyz"):
        atom_types = list(["C"]*self.accessible_points.shape[0]) + list(self.atoms)
        output_pos = np.vstack((self.accessible_points, self.pos))
        ConnollySurface = Atoms(atom_types, output_pos)
        ConnollySurface.write(fname)
        
    def __init__(self, infile):
        martini_radii = pandas.DataFrame()
        self.vdw_radii = pandas.read_csv("Alvarez2013_vdwradii.csv", index_col=0)
        U = mda.Universe(infile)
        self.mol = U.select_atoms("all")
        #mol = read("HCl.pdb")
        self.pos = self.mol.positions
        self.atoms = self.mol.types
        
        self.radius_probe = 1.4
        self.n_sphere_point = 85

        self.sphere_points = self.generate_sphere_points(self.n_sphere_point)
        self.area_per_point = 4.0 * np.pi / len(self.sphere_points) # Scaled in the loop by the vdw_radii
        self.areas = pandas.DataFrame(columns=["area", "atom", "segid", "resname", "resid", "vdw_radius"])

        self.accessible_points = np.ndarray((0, 3))
        
        

if __name__ == "__main__":
    calc = sssSASA("Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe.pdb")
    
    calc.calcSASA()
    
    print(calc.areas)

# =============================================================================
#     for typ in np.unique(calc.areas["atom"]):
#         print(f"Total {typ} area:", calc.areas[calc.areas["atom"] == typ]["area"].sum())
#         print(f"Mean {typ} area:", calc.areas[calc.areas["atom"] == typ]["area"].mean())
# =============================================================================
        
    print(calc.areas["area"].sum())

 
