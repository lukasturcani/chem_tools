query = {}


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


def update(match):
    """
    Adds fingerprints to molecule document.

    Parameters
    ----------
    match : :class:`dict`
        A molecule document from MongoDB.

    Returns
    -------
    :class:`dict`
        A :class:`dict` for updating `match` with
        fingerprints.

    """

    return {

        '$addToSet': {
            'structures.fingerprints'
        }

    }
