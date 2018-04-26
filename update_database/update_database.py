"""
Stores ``stk`` structures in a MongoDB.

"""

import pymongo
import stk
import argparse
from collections import OrderedDict
import numpy as np


def mongo_bb(bb):
    """
    Creates a :class:`dict` to MongoDB with a building block.

    Parameters
    ----------
    bb : :class:`stk.Structunit`
        The building block

    Returns
    -------
    :class:`dict`
        A :class:`dict` which updates MongoDB with `bb` info.

    """

    return {
            '$setOnInsert': {
                              'inchi': bb.inchi,
                              'fg': bb.func_grp.name,
                              'num_fgs': len(bb.functional_group_atoms()),
                              'structure': bb.mdl_mol_block()
                              }
           }


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


def key(molecule):
    """
    Creates a unique key for each `molecule`.

    Paramters
    ---------
    molecule : :class:`stk.Molecule`
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
    parser = argparse.ArgumentParser(description=(
        'Updates a MongoDB with molecular data. \n\n'
        'To identify which molecules to update, the user can provide '
        'either a series of stk Population JSON dump files or a '
        'query file. If stk Population files are provided, then '
        'molecules found in those files are updated. If a query file '
        'is provided, then entries in the MongoDB matching the query '
        'are updated. The query file is a Python file which defines '
        'a function called "query". The function takes two parameters.'
        ' The first one is a molecule and the second is the "key" '
        'function defined in this module. The "query" function must '
        'return a dictionary suitable for use a MongoDB query.'
    ))
    parser.add_argument('mongo_uri', help='URI to the MongoDB server.')
    parser.add_argument('db', help='Name of the database.')
    parser.add_argument('collection', help='Name of the collection.')
    parser.add_argument('update_file',
                        help=('A Python file which defines a '
                              'function called "mongo". The function '
                              'takes 2 parameters. The first is an stk '
                              'Molecule object. The second is the '
                              '"key" function defined in this module. '
                              'The function '
                              'will return a dictionary which is used '
                              'to update the MongoDB.'))

    query_input = parser.add_mutually_exclusive_group()
    query_input.add_argument(
                      '-q',
                      metavar='query_file',
                      help=('A python file which defines a single '
                            'a single variable called "query". This '
                            'holds the query to be used for updating '
                            'MongoDB.'))

    query_input.add_argument(
                        '-p',
                        metavar='population_file',
                        nargs='+',
                        help=('An stk Population JSON dump file. '
                              'The molecules stored in the file will '
                              'be used to update the MongoDB.'))

    args = parser.parse_args()

    client = pymongo.MongoClient(args.mongo_uri)
    col = client[args.db][args.collection]
    bbs = client.bbs.main

    p = stk.Population()
    for pop_file in args.population_file:
        p.add_members(stk.Population.load(pop_file, stk.Molecule.from_dict))

    with open(args.input_file, 'r') as f:
        input_file = {}
        exec(f.read(), {}, input_file)

    for mol in p:
        col.update_one(key(mol), input_file['mongo'](mol, key), True)
        if args.save_bbs:
            for bb in mol.building_blocks:
                bbs.update_one({'inchi': bb.inchi}, mongo_bb(bb), True)


if __name__ == '__main__':
    main()
