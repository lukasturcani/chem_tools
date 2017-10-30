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
c60_files = ['/work/lt912/mmea/databases/targets/C60_OPLS3.pdb',
             '/home/lukas/databases/targets/C60_OPLS3.pdb']
c60_file = next(x for x in c60_files if os.path.exists(x))
settings = {
            'restricted': True,
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

efunc = FunctionData('macromodel', forcefield=16,
                     macromodel_path=mm_path)
ofunc = FunctionData('macromodel_opt', macromodel_path=mm_path,
                     settings=settings, md=md)

fitness_func = FunctionData('cage_c60',
                            target_mol_file=c60_file,
                            efunc=efunc,
                            ofunc=ofunc,
                            n5fold=1, n2fold=1)


class Guard:

    def __init__(self, func_data, dmp):
        self.func = next(x for x in fitness.__dict__.values() if
                         hasattr(x, '__name__') and
                         x.__name__ == func_data.name)
        self.params = func_data.params
        self.dmp = dmp

    def __call__(self, macro_mol):
        try:
            if 'cage_c60' in macro_mol.unscaled_fitness:
                macro_mol.unscaled_fitness.pop('cage_c60')
            if 'cage_c60' in macro_mol.progress_params:
                macro_mol.progress_params.pop('cage_c60')
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
        log_title = 'cage_c60 log - {}.'.format(datetime.now())
        f.write('='*len(log_title) + '\n')
        f.write(log_title + '\n')
        f.write('='*len(log_title) + '\n')
        f.write('Input file: "{}"\n'.format(args.input_file))
        f.write('Output file: "{}"\n'.format(args.output_file))
        f.write('Optimization function: macromodel_opt()\n')
        f.write('settings = {}\n'.format(settings))
        f.write('md = {}\n\n'.format(md))
