from rxn_rebuild.rxn_rebuild import (
    rebuild_rxn,
)
from rxn_rebuild.Args import (
    add_arguments,
    build_args_parser
)
from rr_cache import rrCache
from argparse import (
    ArgumentParser,
    Namespace
)
from logging import (
    Logger
)
from colored import (
    attr as c_attr,
    fg as c_fg,
    bg as c_bg
)


def init(
    parser: ArgumentParser,
    args: Namespace
) -> Logger:
    from brs_utils import create_logger
    from rxn_rebuild._version import __version__

    if args.log.lower() in ['silent', 'quiet'] or args.silent:
        args.log = 'CRITICAL'

    if args.log.lower() in ['silent', 'quiet', 'def_info'] or args.silent:
        disable_rdkit_logging()

    # Create logger
    logger = create_logger(parser.prog, args.log)

    logger.info(
        '{color}{typo}rxn_rebuild {version}{rst}\n'.format(
            prog = logger.name,
            version = __version__,
            color=c_fg('white'),
            typo=c_attr('bold'),
            rst=c_attr('reset')
        )
    )
    logger.debug(args)

    return logger


def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')


def entry_point():
    parser = build_args_parser(
        prog = 'rxn_rebuild',
        description = 'Rebuild full reaction from reaction rule'
    )
    args   = parser.parse_args()

    logger = init(parser, args)

    cache = rrCache(
        db='file',
        attrs=None
    )

    logger.info(
        '{color}{typo}Reaction Rule ID\n   |-{rst} {value}{rst}'.format(
            prog = logger.name,
            value = args.rxn_rule_id,
            color=c_fg('white'),
            typo=c_attr('bold'),
            rst=c_attr('reset')
        )
    )
    logger.info(
        '{color}{typo}Transformation\n   |- to complete:{rst} {value}{rst}'.format(
            prog = logger.name,
            value = args.transfo,
            color=c_fg('white'),
            typo=c_attr('bold'),
            rst=c_attr('reset')
        )
    )
    if args.ori_rxn_id is not None:
        logger.info(
            '{color}{typo}Original (template) reaction ID\n   |-{rst} {value}{rst}'.format(
                prog = logger.name,
                value = args.ori_rxn_id,
                color=c_fg('white'),
                typo=c_attr('bold'),
                rst=c_attr('reset')
            )
        )

    completed_transfos = rebuild_rxn(
              cache = cache,
        rxn_rule_id = args.rxn_rule_id,
            transfo = args.transfo,
        tmpl_rxn_id = args.ori_rxn_id,
             logger = logger
    )

    for tmpl_rxn_id in completed_transfos.keys():
        logger.info(
            '{typo}   |- completed from template reaction {rxn_id}: {rst}{transfo}{rst}'.format(
                prog = logger.name,
                rxn_id = tmpl_rxn_id,
                transfo = completed_transfos[tmpl_rxn_id]['smiles'],
                # color=c_fg('white'),
                typo=c_attr('bold'),
                rst=c_attr('reset')
            )
        )


if __name__ == '__main__':
    entry_point()
