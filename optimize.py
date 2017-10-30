"""
Optimizes molecules held in an ``mtk`` ``.json`` dump file.

This script uses ``MacroModel`` to optimize the molecules.
``MacroModel`` is run using :func:`mtk.macromodel_opt`.

This script is run as a command line program. See::

    $ python optimize.py --help

To run this file a second file with the optimization parameters must
be created. This will be a ``Python`` file which defines two
dictionaries, called `settings` and `md`. These correspond to the
identically named arguments of :func:`mtk.macromodel_opt`. For
example

.. code-block:: python

    # settings_file.py

    settings = {
                'restricted': False,
                'timeout': 0,
                'force_field': 16,
                'max_iter': 2500,
                'gradient': 0.05,
                'md': True
                }

    md = {
          'timeout': 0,
          'force_field': 16,
          'temp': 700,
          'confs': 20,
          'time_step': 1.0,
          'eq_time': 10,
          'sim_time': 200,
          'max_iter': 2500,
          'gradient': 0.05
          }

``optimize.py`` can then be run::

    $ python optimize.py unopt_pop.json settings_file.py opt_pop.json
                         /opt/schrodinger2017-2

"""

from mtk import Population, macromodel_opt
import multiprocessing as mp
from functools import wraps
import os
import argparse
from uuid import uuid4
from os.path import join
from datetime import datetime
import psutil


class Guard:
    """
    Wraps an optimization function.

    This wrapper modifies the function so that its return value is the
    macromolecule it optimizes. It also prevents the function from
    raising errors as this can causes multiprocessing to hang.

    Attributes
    ----------
    mm_path : :class:`str`
        The path to the Schrodinger installation. For example
        ``'/opt/schrodinger2017-2'``.

    settings : :class:`dict`
        A :class:`dict` which holds the arguments for the optimization
        function. Passed to the `settings` argument of the optimization
        function.

    md : :class:`dict`
        A :class:`dict` which holds the arguments for the MD run of
        the optimization function. Passed to the `md` argument of the
        optimization function.

    dmp : :class:`bool`
        If ``True``, a ``.json`` dump file is made for each molecule
        after it is optimized.

    """

    def __init__(self, func, mm_path, settings, md, dmp):
        """
        Initializes a guarded optimization function.

        Parameters
        ----------
        func : :class:`function`
            The optimization function.

        mm_path : :class:`str`
            The path to the Schrodinger installation. For example
            ``'/opt/schrodinger2017-2'``.

        settings : :class:`dict`
            A :class:`dict` which holds the arguments for the
            optimization function. Passed to the `settings` argument of
            the optimization function.

        md : :class:`dict`
            A :class:`dict` which holds the arguments for the MD run of
            the optimization function. Passed to the `md` argument of
            the optimization function.

        dmp : :class:`bool`
            If ``True``, a ``.json`` dump file is made for each
            molecule after it is optimized.

        """

        self.mm_path = mm_path
        self.settings = settings
        self.md = md
        self.dmp = dmp
        wraps(func)(self)

    def __call__(self, macro_mol):
        """
        Runs the optimization function.

        Parameters
        ----------
        macro_mol : :class:`mtk.Molecule`
            The molecule to be optimized.

        Returns
        -------
        :class:`mtk.Molecule`
            The optimized molecule.

        """

        try:
            self.__wrapped__(macro_mol,
                             self.mm_path,
                             self.settings,
                             self.md)

        finally:
            if self.dmp:
                macro_mol.dump(join(self.dmp,
                                    '{}.json'.format(uuid4().int)))
            return macro_mol


if __name__ == '__main__':

    # Set up the command line parser.
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'population_file',
            help=('A .json population dump file. The population'
                  ' holds the molecules which need to be optimized.'))

    parser.add_argument(
            'settings_file',
            help=('A ``.py`` file which defines 2 dictionaries, '
                  '"settings" and "md". These hold the values for the '
                  'equally named arguments in the optimization '
                  'function.'))

    parser.add_argument(
            'output_file',
            help=('The path to a .json dump file. This will hold the'
                  ' population with the optimized molecules.'))

    parser.add_argument(
            'mm_path',
            help='The path to the Schrodinger installation.')

    parser.add_argument(
                '-w', '--write',
                help=('Write .mol files of all optimized molecules'
                      ' into directory WRITE.'))

    parser.add_argument(
            '-d', '--dump',
            help=('After optimizing, a .json file of the molecule'
                  ' is written to the directory DUMP.'))

    parser.add_argument(
            '-n', '--num_cores',
            help='The number of cores to use.',
            type=int,
            default=psutil.cpu_count())

    args = parser.parse_args()

    # Load the data in the settings file.
    with open(args.settings_file, 'r') as f:
        settings_content = {}
        exec(f.read(), settings_content)

    pop = Population.load(args.input_file)
    if args.dump and not os.path.exists(args.dump):
        os.mkdir(args.dump)

    opt_func = Guard(macromodel_opt,
                     args.mm_path,
                     settings_content['settings'],
                     settings_content['md'],
                     args.dump)

    with mp.Pool(args.num_cores) as pool:
        results = pool.map(opt_func, pop)
        if args.write:
            for i, m in enumerate(results):
                m.write(join(args.write, f'{i}.mol'))

    rpop = Population(*results)
    rpop.dump(args.output_file)

    # Write a log of the settings to a file.
    log_file = os.path.splitext(args.output_file)[0] + '.log'
    with open(log_file, 'a') as f:
        log_title = 'Optimization log - {}.'.format(datetime.now())
        f.write('='*len(log_title) + '\n')
        f.write(log_title + '\n')
        f.write('='*len(log_title) + '\n')
        f.write('Input file: "{}"\n'.format(args.input_file))
        f.write('Output file: "{}"\n'.format(args.output_file))
        f.write('Function: macromodel_opt()\n')
        f.write('settings = {}\n'.format(settings_content['settings']))
        f.write('md = {}\n\n'.format(settings_content['md']))
