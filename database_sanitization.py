"""
Holds functions for getting rid of undesirable molecules from databses.

"""

import os
from os.path import join
from glob import iglob
import networkx as nx
import rdkit.Chem.AllChem as rdkit
import itertools as it
import logging
from stk import StructUnit


logger = logging.getLogger(__name__)


def fg_prune(folder, fg, fg_num, ext):
    """
    Deletes molecules without a given functional group from a folder.

    Parameters
    ----------
    folder : :class:`str`
        The full path of the folder holding molecule structure files.
        Any molecules without functional group `fg` `fg_num` amount of
        times are deleted.

    fg : :class:`str`
        The name of the functional group which the desired molecules
        must possess. The name must correspond to one of the name of a
        functional group defined within
        :const:`stk.functional_groups_list`.

    fg_num : :class:`int`
        The number of functional groups of type `fg` which the molecule
        must have in order to not be deleted.

    ext : :class:`str`
        The extension of the molecular structure files in `folder`.
        All other files are skipped.

    Returns
    -------
    None : :class:`NoneType`

    """

    for path in iglob(join(folder, '*.'+ext.lstrip('.'))):
        mol = StructUnit(path, fg)
        # Check that the correct number is present.
        if len(mol.functional_group_atoms()) != fg_num:
            logger.info(f'Deleting {path}.')
            os.remove(path)


def fg_distance_prune(folder, fg, ext):
    """
    Deletes molecules with functional groups seperated by 1 atom.

    Parameters
    ----------
    folder : :class:`str`
        The full path of the folder which holdes the molecule
        structure files. The files are removed from this folder.

    fg : :class:`str`
        The name of the functional group.

    ext : :class:`str`
        The file extension of structure files in `folder`. All other
        files are skipped.

    Returns
    -------
    None : :class:`NoneType`

    """

    for path in iglob(join(folder, '*.'+ext.lstrip('.'))):
        mol = StructUnit(path, fg)
        # Make a mathematical graph of the molecule. Useful for finding
        # the separation between nodes (atoms).
        g = mol.graph()
        # Each bonder atom can act as either a start or end node on a
        # graph. Find the seperations between such nodes. If the
        # separation is 3 this means there is only one node between
        # the start and end nodes. As a result the functional groups
        # are separated by 1 atom and should be deleted.
        for start, end in it.combinations(mol.bonder_ids, 2):
            if nx.shortest_path_length(g, start, end) < 3:
                logger.info(f'Removing {path}.')
                os.remove(path)
                break


def substruct_prune(folder, ext, substruct):
    """
    Deletes molecules which contain the substructure `substruct`.

    Parameters
    ----------
    folder : :class:`str`
        The full path of the folder from which the files are checked
        for the substructure and deleted.

    ext : :class:`str`
        The extension of the structure files in `folder`. All other
        files are skipped.

    substruct : :class:`str`
        The SMARTS string of the substructure, which if present in a
        molecule causes it to be deleted from `folder`.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Create a rdkit molecule of the substructure.
    substruct_mol = rdkit.MolFromSmarts(substruct)

    for path in iglob(join(folder, '*.'+ext.lstrip('.'))):
        mol = StructUnit(path)
        # Check for substruct and delete as appropriate.
        if mol.mol.HasSubstructMatch(substruct_mol):
            logger.info(f'Removing {path}.')
            os.remove(path)
