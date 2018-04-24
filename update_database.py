"""
Stores ``stk`` structures in a MongoDB.

"""

import pymongo
import stk
import argparse
from collections import OrderedDict


def topology_key(topology):
    """

    """

    ...


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
    topology_keys = OrderedDict([(stk.Topology, topology_key)])

    # Take the topology instance of `molecule` and check which
    # function in `topology_keys` to use. The first match is used.
    for cls, key_fn in topology_keys.items():
        if isinstance(molecule.topology, cls):
            top_key_fn = key_fn
            break

    return {
        'class': molecule.__class__.__name__,
        'topology': top_key_fn(molecule.topology),
        'building_blocks': [[bb.__class__.__name__,
                             bb.inchi,
                             bb.func_grp.name] for bb in
                            molecule.building_blocks]
    }



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mongo_uri', help='URI to the MongoDB server.')
    parser.add_argument('db_name', help='Name of the database.')
    parser.add_argument('representation_file',
                        help=('A Python file which defines a single '
                              'function called "mongo". The function '
                              'takes a single parameter, which will '
                              'be an stk Molecule object. The function '
                              'will return a dictionary which can be '
                              'placed into the MongoDB.'))
    parser.add_argument('population_file', nargs='+',
                        help=('An stk Population JSON dump file. '
                              'The molecules stored in the file will '
                              'be saved in the MongoDB.'))

    args = parser.parse_args()
    db = pymongo.MongoClient(args.mongo_uri)[args.db_name]


if __name__ == '__main__':
    main()
