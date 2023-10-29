# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 11:45:50 2023

@author: Alex
"""
import ase
from ase import Atoms
import pandas
import numpy as np
from sklearn.metrics import euclidean_distances
from tqdm import tqdm


class pysasa:
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
    
    def calculate(self, atoms, coordinates, n_sphere_point=24):
        self.atoms = atoms
        self.coordinates = coordinates
        self.accessible_points = np.ndarray((0, 3))
        self.sphere_points = self.generate_sphere_points(n_sphere_point)
        self.area_per_point = 4.0 * np.pi / len(self.sphere_points) # Scaled in the loop by the vdw_radii
        self.areas = pandas.DataFrame(columns=["area", "atom", "vdw_radius"])
        
        mol = Atoms(atoms, coordinates)
        
        #print("Doing hydrogen first!")
        for i in tqdm(np.argsort(mol.numbers)):
        #for i in np.argsort(mol.numbers):
            neighbor_indices = self.find_neighbor_indices(atoms, coordinates, self.radius_probe, i)
            n_neighbor = len(neighbor_indices)
            j_closest_neighbor = 0
            radius = self.radius_probe + self.vdw_radii.at[atoms[i], "vdw_radius"]
            
            n_accessible_point = 0
            for point in self.sphere_points:
                
                is_accessible = True
                # Move sphere point to atomic coord and scale by atomic radius
                test_point = (point*radius) + coordinates[i]
        
                # speed up over np.arange(X, Y) as it starts at a more likely indice meaning less for loop iterations
                # i.e. instead of [0,1,2,3,4,5] it might be [3,4,5,0,1]
                cycled_indices = np.hstack((np.arange(j_closest_neighbor, n_neighbor),
                                            np.arange(0, j_closest_neighbor)))
        
                for j in cycled_indices:
                    coordinates_j = coordinates[neighbor_indices[j]]
                    radius_j = self.radius_probe + self.vdw_radii.at[atoms[j], "vdw_radius"]
                    dist = np.linalg.norm(coordinates_j - test_point)
                    if dist < radius_j:
                        j_closest_neighbor = j
                        is_accessible = False
                        break
                if is_accessible:
                    n_accessible_point += 1
                    self.accessible_points = np.vstack((self.accessible_points, test_point))
                
            area = self.area_per_point*n_accessible_point*radius**2 
            self.areas.loc[i] = [area, atoms[i], self.vdw_radii.at[atoms[i], "vdw_radius"]]
        self.areas = self.areas.sort_index()
        return self.areas["area"].sum()
            
    def writeConnolly(self, fname="ConnollySurface.xyz"):
        atom_types = list(["He"]*self.accessible_points.shape[0]) + list(self.atoms)
        ConnollySurface = Atoms(atom_types,np.vstack((self.accessible_points, self.coordinates)))
        ConnollySurface.write(fname)
        
    def __init__(self, radii_csv):
        self.vdw_radii = pandas.read_csv(radii_csv, index_col=0)
        self.radius_probe = 1.4
        

        