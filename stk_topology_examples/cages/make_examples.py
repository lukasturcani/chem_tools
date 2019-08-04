import stk
import os
import shutil
import argparse
import logging


logger = logging.getLogger(__name__)


def get_args():
    parser = argparse.ArgumentParser(
        description='Create an example cage for each topology.'
    )
    parser.add_argument('output_directory')
    parser.add_argument('--logging_level', default=30, type=int)
    return parser.parse_args()


def main():
    args = get_args()
    logging.basicConfig(level=args.logging_level)

    if os.path.exists(args.output_directory):
        shutil.rmtree(args.output_directory)
    os.mkdir(args.output_directory)

    amine2 = stk.BuildingBlock('NCCN', ['amine'])
    aldehyde2 = stk.BuildingBlock('O=CCCC=O', ['aldehyde'])

    amine3 = stk.BuildingBlock(
        smiles='c1c(N)cc(N)cc1N',
        functional_groups=['amine']
    )
    aldehyde3 = stk.BuildingBlock(
        smiles='c1c(C=O)cc(C=O)cc1C=O',
        functional_groups=['aldehyde']
    )

    amine4 = stk.BuildingBlock(
        smiles='c1c(N)c(N)cc(N)c1N',
        functional_groups=['amine']
    )

    amine5 = stk.BuildingBlock(
        smiles='NC1=C(N)C(N)=C(N)N1N',
        functional_groups=['amine']
    )

    building_blocks = {
        (2, 3): [amine2, aldehyde3],
        (2, 4): [aldehyde2, amine4],
        (2, 5): [aldehyde2, amine5],
        (3, ): [amine3, aldehyde3],
        (3, 4): [aldehyde3, amine4]

    }

    for cls_name, cls in vars(stk.cage).items():
        if (
            not isinstance(cls, type)
            or cls is stk.cage.CageTopology
            or not issubclass(cls, stk.cage.CageTopology)
        ):
            continue

        logger.debug(f'Making {cls_name}.')
        func_groups = set()
        for vertex in cls.vertices:
            func_groups.add(len(vertex.edges))
        func_groups = tuple(sorted(func_groups))

        cage = stk.ConstructedMolecule(
            building_blocks=building_blocks[func_groups],
            topology_graph=cls()
        )

        cage.write(
            os.path.join(args.output_directory, f'{cls_name}.mol')
        )


if __name__ == '__main__':
    main()
