# -*- coding: utf-8 -*-
from ase import Atoms
from ase.io import read
import numpy as np
import pandas


vdw_radii = pandas.read_csv("Alvarez2013_vdwradii.csv", index_col=0)
print(vdw_radii)


# Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe

# =============================================================================
# mol = read("Phe-Phe-Met-Ser-Ile-Arg-Phe-Phe.pdb")
# #ASA
# radius_probe = 1.4
# =============================================================================
# =============================================================================
# nsphere=self.nsphere_var.get()
# atoms = pdbatoms.read_pdb(pdbfile)
# pdbatoms.add_radii(atoms)
# asas = calculate_asa(atoms, radius_probe, nsphere)
# print (asas)
# self.area_var.set(asas)
# =============================================================================
