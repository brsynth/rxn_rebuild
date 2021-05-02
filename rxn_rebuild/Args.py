from argparse  import ArgumentParser
from rxn_rebuild._version import __version__
from typing import(
    Callable,
)


def build_args_parser(
    prog: str,
    description: str = '',
    epilog: str = '',
    m_add_args: Callable = None,
) -> ArgumentParser:

    parser = ArgumentParser(
        prog = prog,
        description = description,
        epilog = epilog
    )

    # Build Parser with rptools common arguments
    parser = add_arguments(parser)

    # Add module specific arguments
    if m_add_args is not None:
        parser = m_add_args(parser)

    return parser


def add_arguments(parser: ArgumentParser) -> ArgumentParser:
    parser.add_argument(
        'rxn_rule_id',
        type=str,
        help='Reaction rule identifier'
    )
    parser.add_argument(
        'transfo',
        type=str,
        help='Transformation to complete with flatten compounds ids (e.g. MNXM181 + MNXM4 = MNXM1 + MNXM1 + MNXM1144) or in SMILES string (e.g. [H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]>>[H]Oc1c([H])c([H])c([H])c([H])c1O[H].O=O.O=O)'
    )
    parser.add_argument(
        '--ori_rxn_id',
        type=str,
        help='Original (template) reaction identifier'
    )
    parser.add_argument(
        '--log', '-l',
        metavar='ARG',
        type=str,
        choices=[
            'debug', 'info', 'warning', 'error', 'critical', 'silent', 'quiet',
            'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'SILENT', 'QUIET'
        ],
        default='def_info',
        help='Adds a console logger for the specified level (default: error)'
    )
    parser.add_argument(
        '--silent', '-s',
        action='store_true',
        default=False,
        help='run %(prog)s silently'
    )
    parser.add_argument(
        '--version', '-v',
        action='version',
        version='%(prog)s {}'.format(__version__),
        help='show the version number and exit'
    )
    return parser
