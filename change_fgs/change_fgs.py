"""
Changes functional group in molecules.

"""

import argparse
import os
import rdkit.Chem.AllChem as rdkit
import multiprocessing as mp
from functools import partial
from glob import glob
from os.path import join
import logging

logger = logging.getLogger(__name__)
# logging.basicConfig(level=logging.DEBUG)

init_funcs = {'.mol': rdkit.MolFromMolFile,
              '.mol2': rdkit.MolFromMol2File}

# When adding new functional groups make sure that the first SMARTS
# in the tuple begins with the atom which has a bond added when the
# fg is being added to another molecule.
fgs = {'amine': ('[N]([H])[H]', '[$([*][N]([H])[H])]'),
       'aldehyde': ('[C](=[O])[H]', '[$([*][C](=[O])[H])]'),
       'carboxylic_acid': ('[C](=[O])-[O][H]',
                           '[$([*][C](=[O])-[O][H])]'),
       'fluorine': ('[F]', '[$([*][F])]'),
       'chlorine': ('[Cl]', '[$([*][Cl])]'),
       'bromine': ('[Br]', '[$([*][Br])]'),
       'iodine': ('[I]', '[$([*][I])]'),
       'alcohol': ('[C]([H])([H])-[O][H]',
                   '[$([*][C]([H])([H])-[O][H])]')}


def flatten(iterable):
    for x in iterable:
        if hasattr(x, '__iter__'):
            yield from flatten(x)
        else:
            yield x


def tag_fg_atoms(mol, fg):
    """
    Adds properties to atoms in `mol` which belong to `fg`.

    Atoms belonging to the functional group are given the property
    'fg'. Atoms not part of the functional group but bonded to it are
    given the property 'attached'.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        The molecule having functional group atoms tagged.

    fg : :class:`tuple` of :class:`str`
        The first string holds the SMARTS identifying atoms to be
        tagged as 'fg'. The second string holds SMARTS identifying
        atoms to be tagged as 'attached'.

    Returns
    -------
    None : :class:`NoneType`

    """

    fg_mol = rdkit.MolFromSmarts(fg[0])
    for idx in flatten(mol.GetSubstructMatches(fg_mol)):
        atom = mol.GetAtomWithIdx(idx)
        atom.SetProp('fg', '1')

    attached_mol = rdkit.MolFromSmarts(fg[1])
    for idx in flatten(mol.GetSubstructMatches(attached_mol)):
        atom = mol.GetAtomWithIdx(idx)
        atom.SetProp('attached', '1')


def remove_fg_atoms(mol):
    """
    Removes all atoms with the tag 'fg'.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        The molecule having functional group atoms tagged.

    Returns
    -------
    :class:`rdkit.Chem.rdchem.Mol`
        The molecule with atoms removed.

    """

    emol = rdkit.EditableMol(mol)
    for atom in reversed(mol.GetAtoms()):
        if atom.HasProp('fg'):
            emol.RemoveAtom(atom.GetIdx())
    return emol.GetMol()


def count_attached(mol):
    """
    Counts the number of atoms labelled 'attached'.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        A molecule to have 'attached' atoms counted.

    Returns
    -------
    :class:`int`
        The number of atoms with the property 'attached' in `mol`.

    """

    return sum(1 for x in mol.GetAtoms() if x.HasProp('attached'))


def bond_fragments(mol):
    """
    Creates bonds beween fragments.

    Creates single bonds between atoms with the property 'attached'.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        A molecule composed of fragments which are to be joined.

    Returns
    -------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        The molecule with all the fragments joined.

    """

    frags = list(rdkit.GetMolFrags(mol))
    main = max(frags, key=len)
    frags.remove(main)
    attached_atoms = [x for x in main if
                      mol.GetAtomWithIdx(x).HasProp('attached')]
    emol = rdkit.EditableMol(mol)
    for main_atom, fg in zip(attached_atoms, frags):
        emol.AddBond(main_atom, fg[0], rdkit.BondType.SINGLE)
    return emol.GetMol()


def add_new_fg(mol, fg):
    """
    The functional group is added to the atoms with the 'attached'
    property.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        The molecule having functional group atoms tagged.

    fg : class:`str`
        The SMARTS of the functional group being added to the molecule.

    Returns
    -------
    :class:`rdkit.Chem.rdchem.Mol`
        The molecule with functional groups added.

    """

    fg = rdkit.MolFromSmarts(fg)
    fg.GetAtomWithIdx(0).SetProp('attached', '1')
    attached = [x for x in mol.GetAtoms() if x.HasProp('attached')]
    for attached_atom in attached:
        mol = rdkit.CombineMols(mol, fg)
    return bond_fragments(mol)


def change_fg(molfile, start, end, fgs):
    """
    Changes the functional group of a molecule.

    Parameters
    ----------
    molfile : :class:`str`
        The path to a .mol or .mol2 molecular structure file.

    start : :class:`str`
        The name of the functional group in the molecule to be changed.

    end : :class:`str`
        The name of the functional group being added to the molecule.

    fgs : :class:`dict`
        Maps names of functional groups to tuples holding SMARTS of
        atoms to be tagged.

    Returns
    -------
    :class:`rdkit.Chem.rdchem.Mol`
        The molecule with the functional group substituted.

    """

    try:
        logger.debug('starting')
        _, ext = os.path.splitext(molfile)
        if ext not in init_funcs:
            raise ValueError(('{} - Not a .mol or'
                              ' .mol2 file.').format(molfile))
        logger.debug('loading')
        mol = rdkit.MolFromMolFile(molfile, removeHs=False)
        logger.debug('tagging')
        tag_fg_atoms(mol, fgs[start])
        logger.debug('removing')
        mol = remove_fg_atoms(mol)
        logger.debug('adding')
        mol = add_new_fg(mol, fgs[end][0])
        # The valence is fucked up unless the Hs get removed and
        # readded.
        mol.RemoveAllConformers()
        mol = rdkit.RemoveHs(mol)
        mol = rdkit.AddHs(mol)
        for atom in mol.GetAtoms():
            atom.UpdatePropertyCache()
        for bond in mol.GetBonds():
            bond.SetStereo(rdkit.BondStereo.STEREONONE)
        rdkit.SanitizeMol(mol)
        rdkit.EmbedMolecule(mol, rdkit.ETKDG())
        logger.debug('done')
        return mol

    except Exception as ex:
        print(ex)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',
                        help=('The molecular structure file of the '
                              'input molecule. Supported file types '
                              'are .mol and .mol2.'))

    parser.add_argument('initial_fg',
                        choices=list(fgs.keys()),
                        help=('The name of the functional group in '
                              'the input molecule.'))

    parser.add_argument('output_file',
                        help=('The molecular structure file to which'
                              ' the substituted molecule is written. '
                              'Must be a .mol file.'))

    parser.add_argument('final_fg',
                        choices=list(fgs.keys()),
                        help=('The name of the functional group in '
                              'the output molecule.'))

    parser.add_argument('-d', '--directory', action='store_true',
                        help=('If set, input_file and output_file are'
                              ' assumed to be directories holding '
                              ' the molecular structure files.'))

    parser.add_argument('-s', '--serial', action='store_true',
                        help=('If set, run serially on directories,'
                              ' rather than in parallel.'))

    args = parser.parse_args()

    if args.directory:
        if not os.path.exists(args.output_file):
            os.mkdir(args.output_file)

        pfunc = partial(change_fg,
                        start=args.initial_fg,
                        end=args.final_fg,
                        fgs=fgs)

        if not args.serial:
            with mp.Pool() as pool:
                new_mols = pool.map(pfunc,
                                    glob(join(args.input_file, '*')))
                logger.debug('applied function')
        else:
            new_mols = []
            for name in glob(join(args.input_file, '*')):
                logger.debug(name)
                new_mols.append(pfunc(name))

        logger.debug('writing')
        for i, mol in enumerate(new_mols):
            rdkit.MolToMolFile(mol,
                               join(args.output_file,
                                    '{}.mol'.format(i)),
                               forceV3000=True)

    else:
        new_mol = change_fg(args.input_file,
                            args.initial_fg,
                            args.final_fg,
                            fgs)
        rdkit.MolToMolFile(new_mol, args.output_file, forceV3000=True)
