"""
Stores ``stk`` structures in a MongoDB.

"""

import pymongo
import stk
import argparse
from collections import OrderedDict
import numpy as np


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
    parser = argparse.ArgumentParser()
    parser.add_argument('mongo_uri', help='URI to the MongoDB server.')
    parser.add_argument('db', help='Name of the database.')
    parser.add_argument('collection', help='Name of the collection.')
    parser.add_argument('input_file',
                        help=('A Python file which defines a single '
                              'function called "mongo". The function '
                              'takes 2 parameters. The first is an stk '
                              'Molecule object. The second is the '
                              '"key" function defined in this module. '
                              'The function '
                              'will return a dictionary which is used '
                              'to update the MongoDB.'))
    parser.add_argument('population_file', nargs='+',
                        help=('An stk Population JSON dump file. '
                              'The molecules stored in the file will '
                              'be saved in the MongoDB.'))

    args = parser.parse_args()

    col = pymongo.MongoClient(args.mongo_uri)[args.db][args.collection]

    p = stk.Population()
    for pop_file in args.population_file:
        p.add_members(stk.Population.load(pop_file, stk.Molecule.from_dict))

    with open(args.input_file, 'r') as f:
        input_file = {}
        exec(f.read(), {}, input_file)

    for mol in p:
        col.update_one(key(mol), input_file['mongo'](mol, key), True)


if __name__ == '__main__':
    main()
