"""
From a two-component cage .mol file, returns:

	1) The SMILES code of the building blocks.
	1.1) The active coordination number of the linkers.
        1.2) The active functional group.
        1.3) The number of times a building block has been used in a cage.

	2) The topology of the cage.

        3) Optional Arguments.
        3.1) Optional (-c): Reform the unrelaxed cage from .mol file.
        3.2) Optional (-b): Reform the unrelaxed building-blocks (bb).

Usage: python reformCage.py <file.mol>

Author: James T. Pegg

"""
import argparse
import stk
from rdkit.Chem import AllChem
import operator

def reform_cage(result, deconstructed_topology):
    """
    Reform the unrelaxed stk cage.

    """
    building_blocks = []
    for k in self:
      bb_mol = AllChem.MolFromSmiles(k)
      bb_Hs = AllChem.AddHs(bb_mol)
      AllChem.EmbedMolecule(bb_Hs)
      if self[k][0] >= 3:
          bb_su = stk.StructUnit3(bb_Hs, [self[k][1]])
          if bb_su not in building_blocks:
              building_blocks.append(bb_su)
      else:
          bb_fg = [self[k][1]]
          bb_su = stk.StructUnit2(bb_Hs, [self[k][1]])
          if bb_su not in building_blocks:
              building_blocks.append(bb_su)

    cage = stk.Cage(building_blocks, deconstructed_topology)
    cage.write(name+'-out.mol')
    cage.dump(name+'-out.json')
    return cage

def write_building_blocks(result):
    """
    Writes the optimised building blocks to .mol file

    """
    wbb_count = 1
    for k in self:
        bb_mol = AllChem.MolFromSmiles(k)
        bb_Hs = AllChem.AddHs(bb_mol)
        AllChem.EmbedMolecule(bb_Hs)
        if self[k][0] >= 3:
            bb_su = stk.StructUnit3(bb_Hs, [self[k][1]])
            bb_su.write(name+'-'+str(wbb_count)+'bb.mol')
        else:
            bb_su = stk.StructUnit2(bb_Hs, [self[k][1]])
            bb_su.write(name+'-'+str(wbb_count)+'bb.mol')
        wbb_count += 1

def topology_calc(coordination_numbers):
    """
    Finds the topology of the system.

    """
    if str(self) == '[(2, 3), (3, 2)]':
        return stk.TwoPlusThree()
    elif str(self) == '[(2, 6), (3, 4)]':
        return stk.FourPlusSix()
    elif str(self) == '[(2, 9), (3, 6)]':
        return stk.SixPlusNine()
    elif str(self) == '[(2, 12), (3, 8)]':
        return stk.EightPlusTwelve()
    elif str(self) == '[(3, 4), (4, 2)]':
        return stk.TwoPlusFour()
    elif str(self) == '[(3, 6), (4, 3)]':
        return stk.ThreePlusSix()
    elif str(self) == '[(3, 12), (4, 6)]':
        return stk.SixPlusTwelve()
    elif str(self) == '[(3, 8)]':
        return stk.FourPlusFour()
    else:
        return 'unknown'

def main():	
    parser = argparse.ArgumentParser(description='Deconstruct cage.mol file')
    parser.add_argument('mol_file')
    parser.add_argument('-c', '--reform_cage',
                        action='store_true',
                        help=('writes unrelaxed cage.'))
    parser.add_argument('-b', '--write_building_blocks',
                        action='store_true',
                        help=('writes unrelaxed building-blocks.'))
    args = parser.parse_args()
    name = str(args.mol_file)[0:-4]

    mol = AllChem.MolFromMolFile(args.mol_file)
    imine = AllChem.MolFromSmiles('CC=NC')
    subs = mol.GetSubstructMatches(imine)

    bonds = []
    for i in subs:
        bonds.append(i[1:3])
    bond_ids = [mol.GetBondBetweenAtoms(x,y).GetIdx() for x,y in bonds]
    new_mol = AllChem.FragmentOnBonds(mol,bond_ids)

    new_smiles = AllChem.MolToSmiles(new_mol, True)
    new_split = new_smiles.split('.')

    building_blocks_dict = {}
    smiles_fragments = []
    uncorrected_smiles_fragments = []
    unknown_functional_group = AllChem.MolFromSmiles('*')

    for i in new_split:
        known_aldehyde = AllChem.ReplaceSubstructs(AllChem.MolFromSmiles(str(i)), AllChem.MolFromSmarts('C=[D1]'), AllChem.MolFromSmarts('C=O'), replaceAll=True)
        unknown_aldehyde = AllChem.ReplaceSubstructs(AllChem.MolFromSmiles(str(i)), AllChem.MolFromSmarts('C=[D1]'), AllChem.MolFromSmarts('C=[D1]'), replaceAll=True)
        smiles_known_aldehyde = AllChem.MolToSmiles(known_aldehyde[0])
        if '*' not in smiles_known_aldehyde:
            smiles_unknown_aldehyde = AllChem.MolToSmiles(unknown_aldehyde[0])
            mol_unknown_aldehyde = AllChem.MolFromSmiles(smiles_unknown_aldehyde)
            functional_group_list = mol_unknown_aldehyde.GetSubstructMatches(unknown_functional_group)
            functional_group_coordination = len(functional_group_list)
            building_blocks_dict[smiles_known_aldehyde] = [functional_group_coordination, 'aldehyde']
            smiles_fragments.append(smiles_known_aldehyde)
            uncorrected_smiles_fragments.append(smiles_unknown_aldehyde)

        else:
            known_amine = AllChem.ReplaceSubstructs(AllChem.MolFromSmiles(str(i)), AllChem.MolFromSmarts('N=[D1]'), AllChem.MolFromSmarts('N'), replaceAll=True)
            unknown_amine = AllChem.ReplaceSubstructs(AllChem.MolFromSmiles(str(i)), AllChem.MolFromSmarts('N=[D1]'), AllChem.MolFromSmarts('[D1]'), replaceAll=True)
            smiles_known_amine = AllChem.MolToSmiles(known_amine[0])
            if '*' not in smiles_known_amine:
                smiles_unknown_amine = AllChem.MolToSmiles(unknown_amine[0])
                mol_unknown_amine = AllChem.MolFromSmiles(smiles_unknown_amine)
                functional_group_list = mol_unknown_amine.GetSubstructMatches(unknown_functional_group)
                functional_group_coordination = len(functional_group_list)
                building_blocks_dict[smiles_known_amine] = [functional_group_coordination, 'amine']
                smiles_fragments.append(smiles_known_amine)
                uncorrected_smiles_fragments.append(smiles_unknown_amine)

    number_of_groups = []
    for i in uncorrected_smiles_fragments:
        a = AllChem.MolFromSmiles(i)
        b = a.GetSubstructMatches(unknown_functional_group)
    number_of_groups.append(len(b))
    number_of_groups_dict = {i:number_of_groups.count(i) for i in number_of_groups}
    coordination_numbers = sorted(number_of_groups_dict.items(), key=operator.itemgetter(0))

    smiles_fragments_number = {i:smiles_fragments.count(i) for i in smiles_fragments}
    result = {key: value + [smiles_fragments_number[key]] for key, value in building_blocks_dict.items()}
    deconstructed_topology = topology_calc(coordination_numbers)

    print(result)
    print(deconstructed_topology)

    if args.reform_cage == True:
        if len(result.keys()) == 2:
        cage = reform_cage(result, deconstructed_topology)
    else:
        print('Error: Three-Component Cage')
    if args.write_building_blocks == True:
        write_building_blocks(result)

if __name__ == '__main__':
  main()
