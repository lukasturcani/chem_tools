"""
Stores ``stk`` structures in a MongoDB.

"""

import pymongo
import stk
import argparse
from collections import OrderedDict
import numpy as np
from textwrap import fill
import multiprocessing as mp


def cage_key(topology):
    """
    Creates a key for :class:`stk.CageTopology`.

    Parameters
    ----------
    topology : :class:`stk.CageTopology`
        A cage topology.

    Returns
    -------
    :class:`dict`
        A key which identifies the topology.

    """

    d = {'class': topology.__class__.__name__}
    d.update({key: list(value) if isinstance(value, np.ndarray) else value
              for key, value in topology.__dict__.items() if
              not key.startswith('_')})
    return d


def topology_key(topology):
    """
    Creates a key for :class:`stk.Topology`.

    Parameters
    ----------
    topology : :class:`stk.Topology`
        A topology.

    Returns
    -------
    :class:`dict`
        A key which identifies the topology.

    """

    d = {'class': topology.__class__.__name__}
    d.update({key: value for key, value in topology.__dict__.items() if
             not key.startswith('_')})
    return d


def struct_unit_key(bb):
    return {'inchi': bb.inchi}


def macromol_key(molecule):
    """
    Creates a unique key for each `molecule`.

    Paramters
    ---------
    molecule : :class:`stk.MacroMolecule`
        The molecule which needs to have a key calculated.

    Returns
    -------
    :class:`dict`
        A key unique to each molecule.

    """

    # A dictionary which maps an stk topology class to a function which
    # creates a key for instances of that class.
    topology_keys = OrderedDict([(stk.CageTopology, cage_key),
                                 (stk.Topology, topology_key)])

    # Take the topology instance of `molecule` and check which
    # function in `topology_keys` to use. The first match is used.
    for cls, key_fn in topology_keys.items():
        if isinstance(molecule.topology, cls):
            top_key_fn = key_fn
            break

    return {
        'class': molecule.__class__.__name__,
        'topology': top_key_fn(molecule.topology),
        'building_blocks': [{'class': bb.__class__.__name__,
                             'inchi': bb.inchi,
                             'fg': bb.func_grp.name} for bb in
                            molecule.building_blocks]
    }


def main():
    description = (
     'Updates a MongoDB.\n\n'
     'There are two major use cases for this script. One where '
     'JSON dump files of stk Populations are used, via the '
     '"-p" option, and one where they are not. In each case the '
     '"input_file" will be slightly different.\n\n'
     'When the "-p" is not used, the input file will define two '
     'things. One will be a variable called "query". This will be a '
     'dictionary which is used to match documents in the database. '
     'The second will be a function called "update". The function '
     'will take one argument. The argument will be a document '
     'matched by the "query" dictionary. The function must return two '
     'things. First, it must return the argument received as input. '
     'Second, it must '
     'return a dictionary, which provides the arguments for '
     'the "update_one" function, with the exception of "filter".\n\n'
     'When the "-p" option is used the input file will only need to '
     'define the "update" function. In this case, the function will '
     'take one parameter. It will be a molecule from the '
     'populations passed via the "-p" option. The function will '
     'return two things. First, the molecule passed as an argument. '
     'Second, a dictionary which provides the arguments for '
     'the "update_one" function, with the exception of "filter".')

    parser = argparse.ArgumentParser(
                description=fill(description, replace_whitespace=False),
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('mongo_uri', help='URI to the MongoDB server.')
    parser.add_argument('db', help='Name of the database.')
    parser.add_argument('collection', help='Name of the collection.')
    parser.add_argument('input_file')
    parser.add_argument('-p',
                        metavar='population_file',
                        dest='population_files',
                        nargs='+',
                        help='An stk Population JSON dump file.')

    args = parser.parse_args()

    client = pymongo.MongoClient(args.mongo_uri)
    col = client[args.db][args.collection]

    with open(args.input_file, 'r') as f:
        exec(f.read(), globals())

    if args.population_files:
        with mp.Pool() as pool:
            p = stk.Population()
            for pop_file in args.population_files:
                p.add_members(stk.Population.load(pop_file,
                                                  stk.Molecule.from_dict))

            update_fn = globals()['update']
            updates = pool.map(update_fn, p)
            for mol, update in updates:
                key = (macromol_key if isinstance(mol, stk.MacroMolecule)
                       else struct_unit_key)
                col.update_one(key(mol), **update)
    else:
        c = col.find(globals()['query'], no_cursor_timeout=True)
        with mp.Pool() as pool:
            update_fn = globals()['update']
            updates = pool.map(update_fn, c)
        for match, update in updates:
            col.update_one({'_id': match['_id']}, **update)
        c.close()


if __name__ == '__main__':
    main()
