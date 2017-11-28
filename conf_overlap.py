import rdkit.Chem.AllChem as rdkit
import numpy as np
import mtk

mtk.CACHE_SETTINGS['ON'] = False

mol1 = mtk.StructUnit('/home/lukas/downloads/c1.mol')
mol2 = mtk.StructUnit('/home/lukas/downloads/c2.mol')

mol1.set_position([0, 0, 0])
mol2.set_position([0, 0, 0])

rdkit.SanitizeMol(mol1.mol)
rdkit.SanitizeMol(mol2.mol)

rmsd, transform = rdkit.GetAlignmentTransform(mol1.mol, mol2.mol)

pmat = mol1.position_matrix()
pmat = np.insert(pmat, 3, 1, 0)
new_pmat = transform @ pmat
new_pmat = np.delete(new_pmat, 3, 0)
mol1.set_position_from_matrix(new_pmat)

rdkit.MolToMolFile(mol1.mol, 'aligned_conf1.mol')
rdkit.MolToMolFile(mol2.mol, 'aligned_conf2.mol')
