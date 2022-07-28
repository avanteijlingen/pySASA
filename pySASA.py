# -*- coding: utf-8 -*-

""" Van der Waals radii in [A] taken from:
A cartography of the van der Waals territories
S. Alvarez, Dalton Trans., 2013, 42, 8617-8636
DOI: 10.1039/C3DT50599E
"""

import MDAnalysis as mda
#from ase.io import read
import numpy as np
import pandas
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from sklearn.metrics import euclidean_distances


#########################
# VARIABLES
radius_probe = 1.4
n_sphere_point = 100
VisualizeSphere = False
#########################
# VALIDATION
#Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe
#http://curie.utmb.edu/getarea.html
#PROBE 1.4 A
# Total  area/energy         =         1366.31
# Number of surface atoms    =           66
# Number of buried atoms     =           12
#########################

def vector(*args):
  """
  Creates a new vector

  Args:
    None: returns zero vector
    v (vector): returns copy of v
    x, y, z: returns (x, y, z)
  """
  n_arg = len(args)
  if n_arg == 0:
    return np.zeros(3, dtype=np.float)
  if n_arg == 1:
    data = args[0]
    if len(data) == 3:
      return np.array(data, copy=True)
    raise TypeError('vector() with 1 argument must have 3 elements')
  if n_arg == 3:
    return np.array(args, dtype=np.float, copy=True)
  else:
    raise TypeError('vector() takes 0, 1 or 3 arguments')
    
def generate_sphere_points(n):
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
      points.append(v3numpy.vector(np.cos(phi)*r, y, np.sin(phi)*r))
   return points

def find_neighbor_indices(atoms, coords, probe, k):
    """
    Returns list of indices of atoms within probe distance to atom k. 
    """
    radius = vdw_radii.at[atoms[k], "vdw_radius"]
    neighbor_indices = []
    d = euclidean_distances(coords[k].reshape(1,3), coords)
    for i in range(d.shape[1]):
        if i == k:
            continue
        radius_i = vdw_radii.at[atoms[i], "vdw_radius"]
        if d[0][i] < radius + radius_i + probe: #+probe twice?
            neighbor_indices.append(i)
    return neighbor_indices

martini_radii = pandas.DataFrame()
vdw_radii = pandas.read_csv("Alvarez2013_vdwradii.csv", index_col=0)
#print(vdw_radii)
#plt.plot(vdw_radii)

# Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe

#mol = read("Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe.pdb")
U = mda.Universe("Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe.pdb")
mol = U.select_atoms("all")
#mol = read("HCl.pdb")
pos = mol.positions
atoms = mol.types

#def calculate_asa(atoms, probe, n_sphere_point):
sphere_points = generate_sphere_points(n_sphere_point)
area_per_point = const = 4.0 * np.pi / len(sphere_points)
areas = pandas.DataFrame(columns=["area", "atom", "segid", "resname", "resid", "vdw_radius"])
result = 0

if VisualizeSphere:
    x = np.vstack(sphere_points)[:,0]
    y = np.vstack(sphere_points)[:,1]
    z = np.vstack(sphere_points)[:,2]
    ax = plt.figure().add_subplot(projection='3d')
    ax.scatter3D(x, y, z, c=z, cmap='Greens');
    plt.show()

for i in range(0, pos.shape[0]):
    neighbor_indices = find_neighbor_indices(atoms, pos, radius_probe, i)
    n_neighbor = len(neighbor_indices)
    j_closest_neighbor = 0
    radius = radius_probe + vdw_radii.at[atoms[i], "vdw_radius"]
    
    n_accessible_point = 0
    for point in sphere_points:
        is_accessible = True
        # Move sphere point to atomic coord and scale by atomic radius
        test_point = (point*radius) + pos[i]

        # speed up over np.arange(X, Y) as it starts at a more likely indice meaning less for loop iterations
        # i.e. instead of [0,1,2,3,4,5] it might be [3,4,5,0,1]
        cycled_indices = np.hstack((np.arange(j_closest_neighbor, n_neighbor),
                                    np.arange(0, j_closest_neighbor)))

        for j in cycled_indices:
            pos_j = pos[neighbor_indices[j]]
            radius_j = radius_probe + vdw_radii.at[atoms[j], "vdw_radius"]
            diff = np.linalg.norm(pos_j - test_point)
            if diff**2 < radius_j**2:
                j_closest_neighbor = j
                is_accessible = False
                break
        if is_accessible:
            n_accessible_point += 1
        
    area = const*n_accessible_point*radius*radius 
    result+=area
    areas.loc[i] = [area, atoms[i], mol[i].segid, mol[i].resname, mol[i].resid, vdw_radii.at[atoms[i], "vdw_radius"]]
print(areas)
    
print(result)

 
