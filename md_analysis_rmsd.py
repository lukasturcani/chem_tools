'''
A simple script that introduces MDAnalysis to 
users for the investigation of molecules.

This example shows how the root-mean-squared displacement
of the C-atoms in a molecule can be calculated.

'''

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

# Specify the reference molecule.
ref = mda.Universe('path-to-file')

# Specify the molecule for comparison.
mobile = mda.Universe('path_to_file')

# Calculates the root-mean-squared displacement
# for the C-skeleton.
result = rmsd(mobile.select_atoms('name C').positions, ref.select_atoms('name C').positions)
print(result)
