"""
Calculates fitness of molecules held in an MMEA .json dump file.

"""

from mmeam import Population, FunctionData
from mmeam.ga import fitness

import multiprocessing as mp
from uuid import uuid4
import os
from os.path import join
import argparse
from datetime import datetime


mm_paths = ['/home/lt912/schrodinger2016-4', '/opt/schrodinger2016-4']
mm_path = next(x for x in mm_paths if os.path.exists(x))


# Options.
pseudoformation_params = {'func': FunctionData('macromodel',
                                               forcefield=16,
                                               macromodel_path=mm_path)
                          }

# A FunctionData instance of the fitness function.
fitness_func = FunctionData(
                        'cage',
                        pseudoformation_params=pseudoformation_params)


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
            'output_file',
            help=('The path to a .json dump file. This will hold the'
                  ' population once the molecules have their fitness'
                  ' values calculation.'))

    parser.add_argument(
                '-w', '--write',
                help=('Write .mol files of all molecules'
                      ' into directory WRITE.'))

    parser.add_argument(
            '-d', '--dump',
            help=('After fitness calculation, a .json file of the '
                  'molecule is written to the directory DUMP.'))

    args = parser.parse_args()
    pop = Population.load(args.input_file)
    if args.dump and not os.path.exists(args.dump):
        os.mkdir(args.dump)

    fit_func = Guard(fitness_func, args.dump)

    with mp.Pool(24) as pool:
        results = pool.map(fit_func, pop)
        if args.write:
            for i, m in enumerate(results):
                m.write(join(args.write, '{}.mol'.format(i)))

    rpop = Population(*results)
    rpop.dump(args.output_file)

    # Write a log of the settings to a file.
    log_file = os.path.splitext(args.output_file)[0] + '.log'
    with open(log_file, 'a') as f:
        log_title = 'cage log - {}.'.format(datetime.now())
        f.write('='*len(log_title) + '\n')
        f.write(log_title + '\n')
        f.write('='*len(log_title) + '\n')
        f.write('Input file: "{}"\n'.format(args.input_file))
        f.write('Output file: "{}"\n\n'.format(args.output_file))
