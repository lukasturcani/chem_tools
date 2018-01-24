"""
Calculates fitness of molecules held in an mtk ``.json`` dump file.

"""

from mtk import Population, Molecule
from mtk.ga import fitness

import multiprocessing as mp
from uuid import uuid4
import os
from os.path import join
import argparse
from datetime import datetime
import psutil


class Guard:

    def __init__(self, func_data, dmp):
        self.func = next(x for x in fitness.__dict__.values() if
                         hasattr(x, '__name__') and
                         x.__name__ == func_data.name)
        self.params = func_data.params
        self.dmp = dmp

    def __call__(self, macro_mol):
        try:
            r = self.func(macro_mol, **self.params)

        except Exception as ex:
            r = None

        finally:
            macro_mol.unscaled_fitness[self.func.__name__] = r
            if self.dmp:
                macro_mol.dump(join(self.dmp,
                                    '{}.json'.format(uuid4().int)))
            return macro_mol


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'input_file',
            help=('A .json population dump file. The population'
                  ' holds the molecules which need to have their '
                  'fitness calcluated.'))

    parser.add_argument(
            'settings_file',
            help=('A ".py" file which defines the fitness function '
                  'as a FunctionData object. This must be saved in '
                  'a "fitness_func" variable.'))

    parser.add_argument(
            'output_file',
            help=('The path to a .json dump file. This will hold the'
                  ' population once the molecules have their fitness'
                  ' values calculated.'))

    parser.add_argument(
                '-w', '--write',
                help=('Write .mol files of all molecules'
                      ' into directory WRITE.'))

    parser.add_argument(
            '-d', '--dump',
            help=('After fitness calculation, a .json file of the '
                  'molecule is written to the directory DUMP.'))

    parser.add_argument(
            '-n', '--num_cores',
            help='The number of cores to use.',
            type=int,
            default=psutil.cpu_count())

    args = parser.parse_args()

    # Load the data in the settings file.
    with open(args.settings_file, 'r') as f:
        settings_namespace = {}
        settings = f.read()
        exec(settings, settings_namespace)

    pop = Population.load(args.input_file, Molecule.fromdict)
    if args.dump and not os.path.exists(args.dump):
        os.mkdir(args.dump)

    fit_func = Guard(settings_namespace['fitness_func'], args.dump)

    with mp.Pool(args.num_cores) as pool:
        results = pool.map(fit_func, pop)
        if args.write:
            for i, m in enumerate(results):
                m.write(join(args.write, f'{i}.mol'))

    rpop = Population(*results)
    rpop.dump(args.output_file)

    # Write a log of the settings to a file.
    log_file = os.path.splitext(args.output_file)[0] + '.log'
    with open(log_file, 'a') as f:
        log_title = 'calc_fitness.py log - {}.'.format(datetime.now())
        f.write('='*len(log_title) + '\n')
        f.write(log_title + '\n')
        f.write('='*len(log_title) + '\n\n')
        f.write(f'Input file: "{args.input_file}"\n')
        f.write(f'Output file: "{args.output_file}"\n')
        f.write(f'Settings file content:\n\n{settings}\n\n')
