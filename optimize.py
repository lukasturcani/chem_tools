"""
Optimizes molecules held in an mtk .json dump file.

For usage run::

    $ python optimize.py --help

"""

from mtkm import Population, macromodel_opt
import multiprocessing as mp
from functools import wraps
import os
import argparse
from uuid import uuid4
from os.path import join
from datetime import datetime

mm_paths = ['/home/lt912/schrodinger2016-4', '/opt/schrodinger2016-4']
mm_path = next(x for x in mm_paths if os.path.exists(x))
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


class Guard:

    def __init__(self, func, mm_path, settings, md, dmp):
        self.mm_path = mm_path
        self.settings = settings
        self.md = md
        self.dmp = dmp
        wraps(func)(self)

    def __call__(self, macro_mol):
        try:
            self.__wrapped__(macro_mol, self.mm_path,
                             self.settings, self.md)

        finally:
            if self.dmp:
                macro_mol.dump(join(self.dmp,
                                    '{}.json'.format(uuid4().int)))

            return macro_mol


if __name__ == '__main__':

    # Set up the command line parser.
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'input_file',
            help=('A .json population dump file. The population'
                  ' holds the molecules which need to be optimized.'))

    parser.add_argument(
            'output_file',
            help=('The path to a .json dump file. This will hold the'
                  ' population with the optimized molecules.'))

    parser.add_argument(
                '-w', '--write',
                help=('Write .mol files of all optimized molecules'
                      ' into directory WRITE.'))

    parser.add_argument(
            '-d', '--dump',
            help=('After optimizing, a .json file of the molecule'
                  ' is written to the directory DUMP.'))

    args = parser.parse_args()

    pop = Population.load(args.input_file)
    if args.dump and not os.path.exists(args.dump):
        os.mkdir(args.dump)

    opt_func = Guard(macromodel_opt, mm_path, settings, md, args.dump)

    with mp.Pool(24) as pool:
        results = pool.map(opt_func, pop)
        if args.write:
            for i, m in enumerate(results):
                m.write(join(args.write, '{}.mol'.format(i)))

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
        f.write('settings = {}\n'.format(settings))
        f.write('md = {}\n\n'.format(md))
