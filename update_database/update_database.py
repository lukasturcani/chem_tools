"""
Stores ``stk`` structures in a MongoDB.

"""

import pymongo
import stk
import argparse
from collections import OrderedDict


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
    d.upate({key: value for key, value in topology.__dict__.items() if
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


def update_databse(db, pop):
    """


    """

    ...


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mongo_uri', help='URI to the MongoDB server.')
    parser.add_argument('db_name', help='Name of the database.')
    parser.add_argument('input_file',
                        help=('A Python file which defines a single '
                              'function called "mongo". The function '
                              'takes two parameters, the first will '
                              'be an stk Molecule object. The second '
                              'parameter is supplied automatically. '
                              'It is the "key" function defined in this '
                              'file. The "mongo" function '
                              'will return a dictionary which is used '
                              'to update the MongoDB.'))
    parser.add_argument('population_file', nargs='+',
                        help=('An stk Population JSON dump file. '
                              'The molecules stored in the file will '
                              'be saved in the MongoDB.'))

    args = parser.parse_args()

    db = pymongo.MongoClient(args.mongo_uri)[args.db_name]

    p = stk.Population()
    for pop_file in args.population_file:
        p.add_members(stk.Population.load(pop_file, stk.Molecule.from_dict))

    update_database(db, p)


if __name__ == '__main__':
    main()
