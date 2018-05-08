def update(molecule):
    """
    Crates a MongoDB document for a molecule.

    This function creates a dictionary which appends a new structure
    to the "structure" field of a molecule in the MongoDB.

    Parameters
    ----------
    molecule : :class:`stk.Molecule`
        The molecule to store in the database.

    Returns
    -------
    :class:`dict`
        Parameters for :func:`update_one`.

    """

    update = {

        '$addToSet': {
                  'tags': {'$each': ['liverpool_refined',
                                     'amines2aldehydes3']},

                  'structures': {
                       'structure': molecule.mdl_mol_block(),
                       'calc_params': {
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
                                  }
                   }
                 }
    return molecule, {'update': update, 'upsert': True}
