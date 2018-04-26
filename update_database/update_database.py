"""
Stores ``stk`` structures in a MongoDB.

"""

import pymongo
import stk
import argparse
from collections import OrderedDict
import numpy as np
from textwrap import fill


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


def su_key(bb):
    return {'inchi': bb.inchi}


def mm_key(molecule):
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
     'will take a single parameter. This will be a document '
     'matched by the "query" dictionary. The function must return '
     'a dictionary, which is used to update the document matched by '
     '"query" and passed as an argument to the function.\n\n'
     'When the "-p" option is used the input file will only need to '
     'define the "update" function. In this case, the function will '
     'take two parameters. The first will be a molecule from the '
     'populations passed via the "-p" option. The second will be '
     'a "key" function defined in this module. The "key" function '
     'will be "su_key" if the molecule is an stk StructUnit object '
     'or "mm_key" if the molecule is an stk MacroMolecule object. '
     'The function will return a dictionary which is used to update '
     'the document corresponding to that molecule in the MongoDB.')

    parser = argparse.ArgumentParser(
                description=fill(description, replace_whitespace=False),
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('mongo_uri', help='URI to the MongoDB server.')
    parser.add_argument('db', help='Name of the database.')
    parser.add_argument('collection', help='Name of the collection.')
    parser.add_argument('input_file')
    parser.add_argument('-p',
                        metavar='population_file',
                        nargs='+',
                        help='An stk Population JSON dump file.')

    args = parser.parse_args()

    client = pymongo.MongoClient(args.mongo_uri)
    col = client[args.db][args.collection]

    with open(args.input_file, 'r') as f:
        input_file = {}
        exec(f.read(), {}, input_file)

    if args.population_file:
        p = stk.Population()
        for pop_file in args.population_file:
            p.add_members(stk.Population.load(pop_file,
                                              stk.Molecule.from_dict))

        for mol in p:
            key = mm_key if isinstance(mol, stk.MacroMolecule) else su_key
            col.update_one(key(mol),
                           input_file['update'](mol, key),
                           True)
    else:
        for match in col.find(input_file['query']):
            col.update_one({'_id': match['_id']},
                           input_file['update'](match),
                           True)


if __name__ == '__main__':
    main()
