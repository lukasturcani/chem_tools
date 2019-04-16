import json
import argparse
import multiprocessing as mp
import stk
import rdkit.Chem.AllChem as rdkit


def load_dict(filename):
    with open(filename, 'r') as f:
        content = json.load(f)
    return content


def write_dict(d, filename):
    with open(filename, 'w') as f:
        json.dump(d, f, indent=4)


def convert_bb_counter(counter):
    return [[convert_bb(bb), count] for bb, count in counter]


def convert_bb(bb):
    new_bb = {'atom_props': {}}
    for key, value in bb.items():
        if key == 'func_grp':
            new_bb['func_groups'] = [value]
        elif key == 'mol_block':
            stk.StructUnit(rdkit.MolFromMolBlock(value, removeHs=False))
            new_bb['conformers'] = [[0, value]]
        else:
            new_bb[key] = value
    return new_bb


def convert(key, value):
    if key == 'mol_block':
        return 'conformers', [[0, value]]
    elif key == 'topology':
        new_value = value.replace('bb_assignments', 'bb_positions')
        return 'topology', new_value
    elif key == 'bb_counter':
        return key, convert_bb_counter(value)
    elif key == 'building_blocks':
        return key, [convert_bb(bb) for bb in value]
    else:
        return key, value


def macromol_func_groups(atom_props):
    func_groups = {}

    for atom, props in atom_props.items():
        if 'fg_id' in props:
            fg_id = props['fg_id']
            if fg_id not in func_groups:
                func_groups[fg_id] = stk.FunctionalGroup(
                                          id_=props['fg_id'],
                                          atom_ids=[],
                                          bonder_ids=[],
                                          deleter_ids=[],
                                          info=props['fg']
                )
            fg = func_groups[fg_id]
            atom_id = int(atom)
            fg.atom_ids.append(atom_id)
            if 'bonder' in props:
                fg.bonder_ids.append(atom_id)

    return repr(list(func_groups.values()))


def convert_macromol(macromol):
    new = {}
    for key, value in macromol.items():
        new_key, new_value = convert(key, value)
        new[new_key] = new_value
    new['func_groups'] = macromol_func_groups(macromol['atom_props'])
    return new


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('old_json')
    parser.add_argument('new_json')
    args = parser.parse_args()

    pop = load_dict(args.old_json)
    with mp.Pool() as pool:
        new_pop = pool.map(convert_macromol, pop)

    write_dict(new_pop, args.new_json)


if __name__ == '__main__':
    main()
