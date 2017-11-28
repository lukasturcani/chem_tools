"""
Checks the angle between the functional groups and centroid of a mol.

"""

import mtk
import numpy as np


def calc_angle(mol_path, fg=None):
    """
    Calculates the angle between the functional groups and centroid.

    Parameters
    ----------
    mol_path : :class:`str`
        The path to a molecular structure file.

    fg : :class:`str`, optional
        The name of the functional group.

    Returns
    -------
    :class:`float`
        The angle between the centroid and functional group atoms in
        degress.

    """

    mol = mtk.StructUnit(mol_path, fg)
    c = mol.centroid()
    a1 = mol.atom_coords(mol.bonder_ids[0]) - c
    a2 = mol.atom_coords(mol.bonder_ids[1]) - c

    rads = np.arccos((a1@a2) / (np.linalg.norm(a1) * np.linalg.norm(a2)))
    return 180*rads / np.pi
