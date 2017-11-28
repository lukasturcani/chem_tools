import rdkit.Chem.AllChem as rdkit
import numpy as np
import mtk
from collections import defaultdict

mtk.CACHE_SETTINGS['ON'] = False


def id_dict(mol):
    d = defaultdict(list)
    for atom in mol.GetAtoms():
        d[atom.GetSymbol()].append(atom.GetIdx())
    return d


def remapped_mol(mol, id_dict):
    emol = rdkit.EditableMol()
    for new_id, old_id in sorted(id_dict.items()):
        emol.AddAtom(mol.GetAtomWithIdx(old_id))

    old_to_new = dict(zip(id_dict.values(), id_dict.keys()))
    for bond in mol.GetBonds():
        id1 = old_to_new[bond.GetBeginAtomIdx()]
        id2 = old_to_new[bond.GetEndAtomIdx()]
        emol.AddBond(id1, id2, bond.GetBondType())

    return emol.GetMol()


def remappings(mol):



def all_rmsds(mol1, mol2):
    sorted1 = sorted_mol(mol1, id_dict1)
    sorted2 = sorted_mol(mol2, id_dict2)


def minimize_rmsd(path1, path2):
    mol1 = mtk.StructUnit(path1)
    mol2 = mtk.StructUnit(path2)

    mol1_ids = id_dict(mol1)
    mol2_ids = id_dict(mol2)

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


path1 = '/home/lukas/downloads/B24_opt_babel1.mol'
path2 = '/home/lukas/downloads/B24_XRD_babel1.mol'
