def bb_fingerprint(molecule, bits, radius):
    """

    """

    ...


def cage_fingerprint(molecule, bits, radius):
    """

    """

    ...


def make_cage_bb_fingerprint(molecule, bits, radius):
    """

    """

    ...


def mongo(molecule, key):
    """

    """

    return {
        '$set': {'key': key(molecule)},

        '$addToSet': {
            'structures.fingerprints'
        }

    }
