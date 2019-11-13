"""
This script uses ``stk`` and ``rdkit`` to find the building-blocks
and topology of a two-component cage molecule held within
a ``.mol`` file.

The script returns:
    1) A dictionary containing:
        Key) The SMILES code of the building blocks.
        Value-1) The coordination number of the linkers.
        Value-2) The active functional group.
        Value-3) Times a building-block has been used in a cage.

    2) The topology of the cage.

    3) Optional Arguments.
        3.1) Optional (-c): Reform the unrelaxed cage from .mol file.
        3.2) Optional (-b): Reform the unrelaxed building-blocks (bb).

This script is run as a command line program. See:

    $ python reform_cage.py <file.mol> --help

Author: James T. Pegg

"""

import argparse
import stk
from rdkit.Chem import AllChem
import operator
import os


def reform_cage(name, building_blocks, deconstructed_topology):
    """
    Reform the unrelaxed cage and write to file.

    Parameters
    ----------
    name : :class:`str`
        The name of the molecule.

    building_blocks_dict : :class:`dict`
        A dictionary containing information on the building blocks.

    deconstructed_topology : :class: ``
        The topology of the molecule.

    """

    building_blocks_list = []
    for k in building_blocks:
        bb_mol = AllChem.MolFromSmiles(k)
        bb_Hs = AllChem.AddHs(bb_mol)
        AllChem.EmbedMolecule(bb_Hs)
        if building_blocks[k][0] >= 3:
            bb_su = stk.StructUnit3(bb_Hs, [building_blocks[k][1]])
            if bb_su not in building_blocks:
                building_blocks_list.append(bb_su)
        else:
            bb_su = stk.StructUnit2(bb_Hs, [building_blocks[k][1]])
            if bb_su not in building_blocks_list:
                building_blocks_list.append(bb_su)

    cage = stk.Cage(building_blocks_list, deconstructed_topology)

    return cage


def write_building_blocks(name, building_blocks):
    """
    Writes the optimised building blocks to ``.mol`` file.

    Parameters
    ----------
    name : :class:`str`
        The name of the molecule.

    building_blocks : :class:`dict`
        A dictionary containing the building blocks.

    Returns
    -------
    None: :class:`NoneType`

    """

    # Create a dictionary to hold the building blocks.
    os.mkdir('BuildingBlocks')
    # Populates the `BuildingBlocks` dictionary.
    wbb_count = 1
    for k in building_blocks:
        bb_mol = AllChem.MolFromSmiles(k)
        bb_Hs = AllChem.AddHs(bb_mol)
        AllChem.EmbedMolecule(bb_Hs)
        if building_blocks[k][0] >= 3:
            bb_su = stk.StructUnit3(bb_Hs, [building_blocks[k][1]])
            bb_su.write('BuildingBlocks/' + name +
                        '-' + str(wbb_count) + 'bb.mol')
        else:
            bb_su = stk.StructUnit2(bb_Hs, [building_blocks[k][1]])
            bb_su.write('BuildingBlocks/' + name +
                        '-' + str(wbb_count) + 'bb.mol')
        wbb_count += 1


def topology_calc(coordination_numbers):
    """
    Finds the topology of the system.

    Parameters
    ----------
    coordination_numbers : :class:`list`
        The number of times a ditopic, tritopic, or quadtopic
        building block occurs.

    Returns
    -------
    The topology of the cage given by `stk`.

    Example
    -------
    The numbers in the list ``[(2, 6), (3, 4)]`` denote:

        (2 = A ditopic building-block.
        6 = The number of times the ditopic building_block occurs.)

        (3 = A tritopic building-block.
        4 = The number of times the tritopic building_block occurs.)

    """

    if coordination_numbers == [(2, 3), (3, 2)]:
        return stk.cage.TwoPlusThree()
    elif coordination_numbers == [(2, 6), (3, 4)]:
        return stk.cage.FourPlusSix()
    elif coordination_numbers == [(2, 9), (3, 6)]:
        return stk.cage.SixPlusNine()
    elif coordination_numbers == [(2, 12), (3, 8)]:
        return stk.cage.EightPlusTwelve()
    elif coordination_numbers == [(2, 4), (4, 2)]:
        return stk.cage.TwoPlusFour()
    elif coordination_numbers == [(2, 6), (4, 3)]:
        return stk.cage.ThreePlusSix()
    elif coordination_numbers == [(2, 12), (4, 6)]:
        return stk.cage.SixPlusTwelve()
    elif coordination_numbers == [(3, 8)]:
        return stk.cage.FourPlusFour()
    elif coordination_numbers == [(3, 8), (4, 6)]:
        return stk.cage.SixPlusEight()
    else:
        return None


def resolve_functional_group(known_smiles,
                             unknown_smiles,
                             functional_group,
                             unknown_functional_group,
                             building_blocks_dict,
                             smiles_fragments,
                             uncorrected_smiles_fragments):
    """
    Edits molecules based on functional group.

    Parameters
    ----------
    known_smiles : :class:`str`
        The smiles of the resolved building block.

    unknown_smiles : :class:`tuple`
        The smiles of the un-resolved building block analogue.

    functional_group : :class:`str`
        The functional group of the building block involved in iminie
        bond formation.

    unknown_functional_group : :class:`rdkit.Chem.rdchem.Mol`
        An unknown `[*]` rdkit mol object needed
        for substructure searches.

    building_blocks_dict : :class:`dict`
        A dictionary containing information on the building blocks.

    smiles_fragments : :class:`list`
        A list containing all the corrected building blocks.

    uncorrected_smiles_fragments : :class:`list`
        A list containing all the uncorrected building blocks.


    Returns
    -------
    None : :class:`NoneType`

    """

    if '*' not in known_smiles:
        unknown_group_smiles = AllChem.MolToSmiles(unknown_smiles[0])
        unknown_group_mol = AllChem.MolFromSmiles(unknown_group_smiles)
        functional_group_list = unknown_group_mol.GetSubstructMatches(
            unknown_functional_group)
        functional_group_coordination = len(functional_group_list)
        building_blocks_dict[known_smiles] = [
            functional_group_coordination, functional_group]
        smiles_fragments.append(known_smiles)
        uncorrected_smiles_fragments.append(unknown_group_smiles)


def main():
    # Set up the command line parser.
    parser = argparse.ArgumentParser()
    parser.add_argument('mol_file')
    parser.add_argument(
            '-c', '--reform_cage',
            action='store_true',
            help=('writes unrelaxed cage.'))

    parser.add_argument(
            '-b', '--write_building_blocks',
            action='store_true',
            help=('writes unrelaxed building-blocks.'))

    args = parser.parse_args()

    # Finds the name of the ``.mol`` file.
    name = str(args.mol_file)[0:-4]

    # Breaks all imine bonds and fragments the molecule.
    mol = AllChem.MolFromMolFile(args.mol_file)
    imine = AllChem.MolFromSmiles('CC=NC')
    subs = mol.GetSubstructMatches(imine)

    bonds = []
    for i in subs:
        bonds.append(i[1:3])
    bond_ids = [mol.GetBondBetweenAtoms(x, y).
                GetIdx() for x, y in bonds]
    new_mol = AllChem.FragmentOnBonds(mol, bond_ids)

    new_smiles = AllChem.MolToSmiles(new_mol, True)
    new_split = new_smiles.split('.')

    # Functionalises the fragmented molecule.
    building_blocks_dict = {}
    smiles_fragments = []
    uncorrected_smiles_fragments = []
    unknown_functional_group = AllChem.MolFromSmiles('*')

    for i in new_split:
        known_aldehyde = AllChem.ReplaceSubstructs(
            mol=AllChem.MolFromSmiles(str(i)),
            query=AllChem.MolFromSmarts('C=[D1]'),
            replacement=AllChem.MolFromSmarts('C=O'),
            replaceAll=True
        )
        unknown_aldehyde = AllChem.ReplaceSubstructs(
            mol=AllChem.MolFromSmiles(str(i)),
            query=AllChem.MolFromSmarts('C=[D1]'),
            replacement=AllChem.MolFromSmarts('C=[D1]'),
            replaceAll=True
        )
        smiles_known_aldehyde = AllChem.MolToSmiles(known_aldehyde[0])
        resolve_functional_group(
            known_smiles=smiles_known_aldehyde,
            unknown_smiles=unknown_aldehyde,
            functional_group='aldehyde',
            unknown_functional_group=unknown_functional_group,
            building_blocks_dict=building_blocks_dict,
            smiles_fragments=smiles_fragments,
            uncorrected_smiles_fragments=uncorrected_smiles_fragments
        )

        known_amine = AllChem.ReplaceSubstructs(
            mol=AllChem.MolFromSmiles(str(i)),
            query=AllChem.MolFromSmarts('N=[D1]'),
            replacement=AllChem.MolFromSmarts('N'),
            replaceAll=True
        )
        unknown_amine = AllChem.ReplaceSubstructs(
            mol=AllChem.MolFromSmiles(str(i)),
            query=AllChem.MolFromSmarts('N=[D1]'),
            replacement=AllChem.MolFromSmarts('[D1]'),
            replaceAll=True
        )
        smiles_known_amine = AllChem.MolToSmiles(known_amine[0])
        resolve_functional_group(
            known_smiles=smiles_known_amine,
            unknown_smiles=unknown_amine,
            functional_group='amine',
            unknown_functional_group=unknown_functional_group,
            building_blocks_dict=building_blocks_dict,
            smiles_fragments=smiles_fragments,
            uncorrected_smiles_fragments=uncorrected_smiles_fragments
        )

    # Collect information.
    number_of_groups = []
    for i in uncorrected_smiles_fragments:
        a = AllChem.MolFromSmiles(i)
        b = a.GetSubstructMatches(unknown_functional_group)
        number_of_groups.append(len(b))
        number_of_groups_dict = {
            i: number_of_groups.count(i) for i in number_of_groups}
        coordination_numbers = sorted(
            number_of_groups_dict.items(),
            key=operator.itemgetter(0))

    smiles_fragments_number = {
        i: smiles_fragments.count(i) for i in smiles_fragments}
    building_blocks = {key: value + [smiles_fragments_number[key]]
                       for key, value in building_blocks_dict.items()}
    deconstructed_topology = topology_calc(coordination_numbers)

    # Prints the output.
    print(building_blocks)
    print(deconstructed_topology)

    if args.reform_cage:
        if len(building_blocks.keys()) == 2:
            cage = reform_cage(name, building_blocks, deconstructed_topology)
            cage.write(name + '-out.mol')
            cage.dump(name + '-out.json')
        else:
            raise Exception('Error: Three-Component Cage')

    if args.write_building_blocks:
        write_building_blocks(name, building_blocks)


if __name__ == '__main__':
    main()

