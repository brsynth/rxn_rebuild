from rxn_rebuild.rxn_rebuild import (
    rebuild_rxn,
)
from rxn_rebuild.Args import build_args_parser
from rr_cache import rrCache
from argparse import (
    ArgumentParser,
    Namespace
)
from logging import (
    Logger,
    basicConfig as log_basicConfig,
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
    if args.log_file != '':
        log_basicConfig(filename=args.log_file, encoding='utf-8')

    cache = rrCache(
        attrs=['rr_reactions', 'template_reactions', 'cid_strc']
        # logger=logger
    )

    msg_rr = '{color}{typo}Reaction Rule\n   |- ID:{rst} {rr_id}'
    if args.tmpl_rxn_id is not None:
        msg_rr += '\n{color}{typo}   |- template reaction:{rst} {tmpl_rxn_id}'
    logger.info(
        msg_rr.format(
            rr_id = args.rxn_rule_id,
            tmpl_rxn_id = args.tmpl_rxn_id,
            color=c_fg('white'),
            typo=c_attr('bold'),
            rst=c_attr('reset')
        )
    )
    logger.info(
        '{color}{typo}Transformation\n   |- to complete:{rst} {value}'.format(
            value = args.transfo,
            color=c_fg('white'),
            typo=c_attr('bold'),
            rst=c_attr('reset')
        )
    )

    completed_transfos = rebuild_rxn(
              cache = cache,
        rxn_rule_id = args.rxn_rule_id,
            transfo = args.transfo,
          direction = args.direction,
        tmpl_rxn_id = args.tmpl_rxn_id,
             logger = logger
    )

    for tmpl_rxn_id in completed_transfos.keys():
        if 'full_transfo' in completed_transfos[tmpl_rxn_id]:
            logger.info(
                '{typo}   |- completed from template reaction {rxn_id}: {rst}{transfo}'.format(
                    rxn_id=tmpl_rxn_id,
                    transfo=completed_transfos[tmpl_rxn_id]['full_transfo'],
                    typo=c_attr('bold'),
                    rst=c_attr('reset')
                )
            )
        else:
            logger.info(
                '{typo}   |- completed from template reaction {rxn_id}: {rst}'.format(
                    rxn_id=tmpl_rxn_id,
                    typo=c_attr('bold'),
                    rst=c_attr('reset')
                )
            )
            logger.info(
                '{typo}{color}         Unknown structure for some compounds{rst}'.format(
                    typo=c_attr('bold'),
                    color=c_fg('white'),
                    rst=c_attr('reset')
                )
            )
            for side in ['left', 'right']:
                logger.info(
                    '{rst}{typo}            |- {side}: {rst}{compounds}{rst}'.format(
                        side=side.upper(),
                        compounds=' '.join(list(completed_transfos[tmpl_rxn_id]['added_cmpds'][side+'_nostruct'].keys())),
                        typo=c_attr('bold'),
                        rst=c_attr('reset')
                    )
                )


if __name__ == '__main__':
    entry_point()
