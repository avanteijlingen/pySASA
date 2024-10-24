# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 11:46:12 2023

@author: Alex
"""

import pysasa
from ase import Atoms
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import sys, pandas
from sklearn.metrics import euclidean_distances

# =============================================================================
# Validation
# Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe
# http://curie.utmb.edu/getarea.html
#  Probe radius  :  1.400
#   Residue      Total   Apolar  Backbone Sidechain Ratio(%) In/Out
#  PHE      1   223.07   189.49    41.18   181.90    100.0     o
#  PHE      2   174.90   172.93     4.54   170.37     94.6     o
#  MET      3   170.00   168.43     6.09   163.92    100.0     o
#  SER      4    64.76    30.58     9.09    55.67     71.9     o
#  ILE      5   100.21    98.04     8.00    92.21     62.6     o
#  ARG      6   200.13    91.81    22.59   177.53     90.8     o
#  PHE      7   199.90   181.80    26.61   173.30     96.2     o
#  PHE      8   233.33   178.96    41.39   191.94    100.0     o
#  ---------------------------------------------
#  POLAR  area/energy         =          254.28
#  APOLAR area/energy         =         1112.03
#  UNKNOW area/energy         =            0.00
#  ---------------------------------------------
#  Total  area/energy         =         1366.31
#  ---------------------------------------------
#  Number of surface atoms    =           66
#  Number of buried atoms     =           12
#  Number of atoms with ASP=0 =           76
# =============================================================================


#mol = read("ExampleData/Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe.pdb")
mol = read("ExampleData/Vignesh_Rendering_100ps/II_wrapped.pdb")
#mol = read("ExampleData/CF3CH3.xyz")
#mol = read("C:/Users/Alex/Documents/GitHub/ANI-based-carbene-pKa-prediction/qmspin_extra/Training18.xyz")

radius_probe = 3

calc = pysasa.pysasa(radii_csv="ExampleData/Alvarez2013_vdwradii.csv", radius_probe = radius_probe)

sasa = calc.calculate(mol.get_chemical_symbols(), mol.positions, n_sphere_point=50)
print("\n")
print("SASA:", sasa)

print(calc.areas)

calc.writeConnolly(f"ConnollySurface_II_{radius_probe}.xyz")

sys.exit()


# =============================================================================
# calc = pysasa.pysasa(radii_csv="ExampleData/orca_vdwradii.csv")
# 
# U = mda.Universe("ExampleData/Carbenes/1_deprot.xyz")
# protein = U.select_atoms("all")
# 
# d = euclidean_distances(protein.positions.astype(np.float16), protein.positions.astype(np.float16))
# 
# 
# #mol = read("ExampleData/Alvaro/PartEq.gro")
# mol = Atoms([x[0] for x in protein.names], protein.positions)
# names = np.array(mol.get_chemical_symbols(), dtype="<U2")
# # Extract just the protein
# mol = mol[np.where((names == "B") | (names == "S"))[0]]
# 
# 
# 
# calc.calculate(mol.get_chemical_symbols(), mol.positions, n_sphere_point=24)
# =============================================================================


calc = pysasa.pysasa(radii_csv="ExampleData/orca_vdwradii.csv")
for state in ["prot", "deprot"]:
    SASA = []
    number_of_atoms = []
    for i in range(1, 12):
        if state == "deprot":
            mol = read(f"ExampleData/Carbenes/{i}.xyz")
        else:
            mol = read(f"ExampleData/Carbenes/{i}+.xyz")
        
        sasa = calc.calculate(mol.get_chemical_symbols(), mol.positions, n_sphere_point=24)
        print()
        print(i, sasa)
        print(calc.areas)
        SASA.append(sasa)
        number_of_atoms.append(mol.positions.shape[0])
        plt.text(number_of_atoms[-1], SASA[-1], str(i))
    plt.scatter(number_of_atoms, SASA, label="state")
plt.ylabel("Surface area (DFT OPT)")
plt.xlabel("Number of Atoms")
