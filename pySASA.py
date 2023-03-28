# -*- coding: utf-8 -*-

""" Van der Waals radii in [A] taken from:
A cartography of the van der Waals territories
S. Alvarez, Dalton Trans., 2013, 42, 8617-8636
DOI: 10.1039/C3DT50599E
"""

import MDAnalysis as mda
from ase import Atoms
import numpy as np
import pandas, time, sys, os
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from sklearn.metrics import euclidean_distances
from tqdm import tqdm

#########################
# VALIDATION
#Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe
#http://curie.utmb.edu/getarea.html
#PROBE 1.4 A
# Total  area/energy         =         1366.31
# Number of surface atoms    =           66
# Number of buried atoms     =           12
#########################


        
class SASACRUNCH:
    def Fibb(self, n):
        goldenRatio = (1 + 5**0.5)/2 #\phi = golden ratio
        i = np.arange(0, n)
        theta = 2 *np.pi * i / goldenRatio
        phi = np.arccos(1 - 2*(i+0.5)/n)
        x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
        #print("theta", theta, theta.shape[0])
        #print("phi", phi, phi.shape[0])
        points = np.array((x,y,z)).T
        return points
    def generate_sphere_points(self, n):
       """
       Returns list of coordinates on a sphere using the Golden-
       Section Spiral algorithm.
       """
       points = np.ndarray((n, 3))
       inc = np.pi * (3 - np.sqrt(5))
       offset = 2 / float(n)
       for k in range(int(n)):
          y = k * offset - 1 + (offset / 2)
          r = np.sqrt(1 - y*y)
          phi = k * inc
          point = np.array((np.cos(phi)*r, y, np.sin(phi)*r), dtype=np.float64, copy=True)
          points[k] = point
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
        
    def calcSASA(self, Frame=-1):
        self.U.trajectory[Frame]
        st = time.time()
        for i in tqdm(range(0, self.pos.shape[0])):
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
                    dist = np.linalg.norm(pos_j - test_point)
                    if dist < radius_j:
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
    

    def double_cubic_lattice_method(self):
        self.atomtic_distances = euclidean_distances(self.pos, self.pos)
        self.radii = self.radius_probe + self.vdw_radii.loc[self.atoms]["vdw_radius"].values
        
        #There must be a faster way to do the following within numpy
        self.surfaces = np.ndarray((self.radii.shape[0], self.n_sphere_point, 3))
        for i in range(self.radii.shape[0]):
            self.surfaces[i] = (self.sphere_points * self.radii[i]) + self.pos[i]
        self.surfaces = self.surfaces.reshape(-1, 3)
        
        #Remove points that are inaccessible
        self.surface_distances = euclidean_distances(self.surfaces, self.surfaces)
        np.fill_diagonal(self.surface_distances, 100)
        
        self.accessible = np.ones((self.surfaces.shape[0]), dtype=np.bool_)
        
        #can get the table of sums by multiplying, deviding
        self.cutoffs = np.add.outer(self.radii, self.radii)
        self.cutoffs += self.radius_probe
        
        for i in range(self.surface_distances.shape[0]):
            self.cuff = self.cutoffs[np.floor(i/self.pos.shape[0]).astype(np.int64)]
            break
        
        
# =============================================================================
#         for point in self.sphere_points:
#             is_accessible = True
#             # Move sphere point to atomic coord and scale by atomic radius
#             test_point = (point*radius) + self.pos[i]
#             print(test_point)
#             break
# =============================================================================
        

class sssSASA(SASACRUNCH): 
    def writeConnolly(self, fname="ConnollySurface.xyz"):
        atom_types = list(["C"]*self.accessible_points.shape[0]) + list(self.atoms)
        output_pos = np.vstack((self.accessible_points, self.pos))
        ConnollySurface = Atoms(atom_types, output_pos)
        ConnollySurface.write(fname)
        
    def VisSphere(self):
        x = np.vstack(self.sphere_points)[:,0]
        y = np.vstack(self.sphere_points)[:,1]
        z = np.vstack(self.sphere_points)[:,2]
        ax = plt.figure().add_subplot(projection='3d')
        #ax.scatter3D(x, y, z, c=z, cmap='Greens')
        ax.scatter3D(*self.sphere_points.T)
        plt.show()
        
    def __init__(self, radii_file,  infiles, n_sphere_point = 150, sel = "not resname H2O and not resname W and not resname CL and not resname NA"):
        #martini_radii = pandas.DataFrame()
        self.vdw_radii = pandas.read_csv(radii_file, index_col=0)
        self.U = mda.Universe(*infiles)
        self.mol = self.U.select_atoms(sel)
        self.pos = self.mol.positions
        self.atoms = self.mol.types
        
        self.radius_probe = 1.4
        self.n_sphere_point = n_sphere_point

        self.sphere_points = self.generate_sphere_points(self.n_sphere_point)
        self.area_per_point = 4.0 * np.pi / len(self.sphere_points) # Scaled in the loop by the vdw_radii
        self.areas = pandas.DataFrame(columns=["area", "atom", "segid", "resname", "resid", "vdw_radius"])

        self.accessible_points = np.ndarray((0, 3))
        self.Fibb_points = self.Fibb(self.n_sphere_point)

if __name__ == "__main__":
    
    jobs = ["FF_2.1/FF_min"]
    #jobs = ["FF_2.1/FF_eq_centred"]
    
    for job in jobs:
        oname = f"{job}.csv"
# =============================================================================
#         if os.path.exists(f"{oname}.csv"):
#             print("Exists:", f"{oname}.csv")
#             continue
# =============================================================================
        print(oname)
        
        radii_file = "Martini_vdwradii.csv"
        #radii_file = "Alvarez2013_vdwradii.csv"
        
        #calc = sssSASA(infiles = [f"{job}.psf", f"{job}.pdb"], n_sphere_point=10)
        calc = sssSASA(radii_file, infiles = [f"{job}.gro"], n_sphere_point=24)
        
        a="""
        for typ in np.unique(calc.U.atoms.types):
            if typ[0] == "S":
                print(f"{typ},4.3")
            else:
                print(f"{typ},4.7")
        #"""
        
        frame = 0

        calc.calcSASA(Frame = frame)
        print()
        calc.areas.to_csv(oname)
        
        
        print(calc.areas["area"].sum())
    
    
        calc.writeConnolly()
    
    #print(calc.areas)

# =============================================================================
#     for typ in np.unique(calc.areas["atom"]):
#         print(f"Total {typ} area:", calc.areas[calc.areas["atom"] == typ]["area"].sum())
#         print(f"Mean {typ} area:", calc.areas[calc.areas["atom"] == typ]["area"].mean())
# =============================================================================

