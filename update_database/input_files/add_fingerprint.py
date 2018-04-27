import rdkit.Chem.AllChem as rdkit
import itertools as it

query = {}


def bb_fingerprint(mols, bits, radius):
    """

    """

    bbs = mols[1:]
    full_fp = []
    for mol in bbs:
        fp = rdkit.GetMorganFingerprintAsBitVect(mol, radius, bits)
        full_fp.extend(list(fp))
    return full_fp


def cage_fingerprint(mols, bits, radius):
    """

    """

    return list(rdkit.GetMorganFingerprintAsBitVect(mols[0], radius, bits))


def cage_bb_fingerprint(mols, bits, radius):
    """

    """

    full_fp = []
    for mol in mols:
        fp = rdkit.GetMorganFingerprintAsBitVect(mol, radius, bits)
        full_fp.extend(list(fp))
    return full_fp


def update(match, client):
    """
    Adds fingerprints to molecule document.

    Parameters
    ----------
    match : :class:`dict`
        A molecule document from MongoDB.

    client : :class:`MongoClient`
        The database client.

    Returns
    -------
    :class:`dict`
        Parameters for :func:`update_one`.

    """

    calc_params = {
                    'software': 'schrodinger2017-4',
                    'max_iter': 5000,
                    'md': {
                           'confs': 50,
                           'temp': 700,
                           'sim_time': 2000,
                           'time_step': 1.0,
                           'force_field': 16,
                           'max_iter': 2500,
                           'eq_time': 100,
                           'gradient': 0.05,
                           'timeout': 0
                    },
                    'force_field': 16,
                    'restricted': 'both',
                    'timeout': 0,
                    'gradient': 0.05
    }

    cage_block = next(structure['structure'] for structure
                      in match['structures'] if
                      structure['calc_params'] == calc_params)

    bb1_inchi = match['key']['building_blocks'][0]['inchi']
    bb1_block = client.bbs.main.find_one({'inchi': bb1_inchi})['structure']
    bb2_inchi = match['key']['building_blocks'][1]['inchi']
    bb2_block = client.bbs.main.find_one({'inchi': bb2_inchi})['structure']

    cage = rdkit.MolFromMolBlock(cage_block, sanitize=False)
    bb1 = rdkit.MolFromMolBlock(bb1_block, sanitize=False)
    bb2 = rdkit.MolFromMolBlock(bb2_block, sanitize=False)
    mols = (cage, bb1, bb2)

    for mol in mols:
        rdkit.GetSSSR(mol)
        mol.UpdatePropertyCache(strict=False)

    radii = [1, 2, 4, 8, 16]
    bit_sizes = [256, 512, 1024, 2048]
    featurizations = [('bb', bb_fingerprint),
                      ('cage', cage_fingerprint),
                      ('cage+bb', cage_bb_fingerprint)]

    fingerprints = [{'type': ['rdkit', 'morgan', name],
                     'bits': bits,
                     'radius': radius,
                     'fp': fn(mols, bits=bits, radius=radius)} for
                    radius, bits, (name, fn) in
                    it.product(radii, bit_sizes, featurizations)]

    update = {

        '$addToSet': {
            'structures.$[s].fingerprints': {'$each': fingerprints}
         }

    }

    array_filters = [{'s.calc_params': calc_params}]

    return {'update': update,
            'upsert': True,
            'array_filters': array_filters}
