from logging import (
    Logger,
    getLogger,
)
from typing import (
    List,
    Dict,
    Tuple
)
from collections import Counter
from json import dumps
from copy import deepcopy
from rr_cache import rrCache
from chemlite import Reaction
from .Args import DEFAULTS


def rebuild_rxn(
    rxn_rule_id: str,
    transfo: str,
    tmpl_rxn_id: str = None,
    cache: 'rrCache' = None,
    cmpds_to_ignore: List[str] = [],
    cspace: str = DEFAULTS['cspace'],
    cspace_type: str = 'rr2026',
    logger: Logger = getLogger(__name__)
) -> str:

    logger.debug(f'rxn_rule_id: {rxn_rule_id}')
    logger.debug(f'transfo: {transfo}')
    logger.debug(f'tmpl_rxn_id: {tmpl_rxn_id}')
    logger.debug(f'cmpds_to_ignore: {cmpds_to_ignore}')
    logger.debug(f'cspace: {cspace}')
    logger.debug(f'cspace_type: {cspace_type}')

    ## INPUT TRANSFORMATION
    trans_input = Reaction.parse(transfo, logger)

    ## LOAD CACHE
    if cache is None:
        # cache = rrCache(
        #     attrs=['rr_reactions', 'template_reactions', 'cid_strc']
        #     # logger=logger
        # )
        cache = rrCache(
            cspace=cspace,
            interactive=False,
            logger=logger
        )

    ## COMPLETE TRANSFORMATION
    try:
        legacy = cspace_type=='legacy'
        completed_transfos = {}
        if not tmpl_rxn_id is None:
            completed_transfos[tmpl_rxn_id] = complete_transfo(
                trans_input=trans_input,
                rxn_rule=cache.get('rr_reactions')[rxn_rule_id][tmpl_rxn_id],
                tmpl_rxn=cache.get('template_reactions')[tmpl_rxn_id],
                tmpl_rxn_id=tmpl_rxn_id,
                compounds=cache.get('cid_strc'),
                cmpds_to_ignore=cmpds_to_ignore,
                legacy=legacy,
                logger=logger
            )
        else:  # One completed transformation per template reaction
            for tpl_rxn_id, rxn_rule in cache.get('rr_reactions')[rxn_rule_id].items():
                completed_transfos[tpl_rxn_id] = complete_transfo(
                    trans_input=trans_input,
                    rxn_rule=rxn_rule,
                    tmpl_rxn=cache.get('template_reactions')[tpl_rxn_id],
                    tmpl_rxn_id=tpl_rxn_id,
                    compounds=cache.get('cid_strc'),
                    cmpds_to_ignore=cmpds_to_ignore,
                    legacy=legacy,
                    logger=logger
                )
    except KeyError as e:
        logger.error(f'   |- KeyError: {str(e)}')
        logger.error('      + The reaction rule is not known in the cache. Are you sure you provided the right data-type, e.g. mnx3.1, mnx4.4...?')
        return {}

    return completed_transfos


def complete_transfo(
    trans_input: Dict,
    rxn_rule: Dict,
    tmpl_rxn: Dict,
    tmpl_rxn_id: str,
    compounds: Dict,
    cmpds_to_ignore: List[str] = [],
    legacy: bool = False,
    logger: Logger = getLogger(__name__)
) -> Dict:

    logger.debug(f'TRANS_INPUT: {dumps(trans_input, indent=4)}')
    logger.debug(f'REACTION RULE: {dumps(rxn_rule, indent=4)}')
    logger.debug(f'TEMPLATE REACTION ({tmpl_rxn_id}): {dumps(tmpl_rxn, indent=4)}')
    # logger.debug(f'COMPOUNDS: {dumps(compounds, indent=4)}')
    logger.debug(f'CMPDS TO IGNORE: {cmpds_to_ignore}')

    ## CHECK 1/2
    # Check if the number of structures in the right part of SMILES of transformation to complete
    # is equal to the number of products of the template reaction used to build the reaction rule.
    # Just in right part since rules are always mono-substrate
    check_compounds_number(
        'INPUT TRANSFORMATION [right]',
        trans_input['right'],
        'REACTION RULE [right]',
        rxn_rule['right'],
        logger=logger
    )

    if legacy:
        logger.debug(f'Entering in legacy rules mode')

        ## TEMPLATE REACTION
        if rxn_rule['rel_direction'] == -1:
            _tmpl_rxn = {
                'right': tmpl_rxn['left'],
                'left': tmpl_rxn['right']
            }
        else:
            _tmpl_rxn = {
                'right': tmpl_rxn['right'],
                'left': tmpl_rxn['left']
            }

        ## ADD MISSING COMPOUNDS TO THE FINAL TRANSFORMATION
        # to replace by LEFT_EXCLUDEED_IDS and RIGHT_EXCLUDED_IDS from templates.tsv
        missing_compounds = detect_missing_compounds(
            _tmpl_rxn,
            rxn_rule,
            compounds,
            cmpds_to_ignore=cmpds_to_ignore,
            logger=logger
        )
    else:
        logger.debug(f'Entering in new rules mode')
        missing_compounds = {
            'left': dict((Counter(rxn_rule['left_excluded']))),
            'right': dict((Counter(rxn_rule['right_excluded'])))
        }
    logger.debug('MISSING COMPOUNDS: '+str(dumps(missing_compounds, indent=4)))

    ## BUILD FINAL TRANSFORMATION
    compl_transfo = build_final_transfo(
        trans_input=trans_input,
        missing_compounds=missing_compounds,
        logger=logger
    )

    ## CHECK 2/2
    # Check if the number of compounds in both right and left sides of SMILES of the completed transformation
    # is equal to the ones of the template reaction
    for side in Reaction.get_SIDES():
        # Adjust template reaction side according to rule relative direction
        # _side = side
        _side = ('left' if side == 'right' else 'right') if rxn_rule.get('rel_direction') == -1 else side
        _tmpl_rxn_side = tmpl_rxn[_side]
        check = check_compounds_number(
            f'COMPLETED TRANSFORMATION ({rxn_rule["rule_id"]}) [{side}]',
            compl_transfo[side],
            f'TEMPLATE REACTION ({tmpl_rxn_id}) [{_side}]',
            _tmpl_rxn_side,
            logger=logger
        )
        if not check:
            exit()

    return {
        'full_transfo': compl_transfo,
        'added_cmpds': missing_compounds,
        'sep_side': trans_input['sep_side'],
        'sep_cmpd': trans_input['sep_cmpd']
    }


def check_compounds_number(
    rxn_name_1: str,
    rxn_1_side: Dict,
    rxn_name_2: str,
    rxn_2_side: Dict,
    logger: Logger = getLogger(__name__)
) -> bool:
    logger.debug(f'Checking number of compounds between {rxn_name_1} and {rxn_name_2}...')
    logger.debug(f'{rxn_name_1}: {rxn_1_side}')
    logger.debug(f'{rxn_name_2}: {rxn_2_side}')
    # Check if the number of structures in the part of SMILES of rxn_1
    # is equal to the number of products of rxn_2.
    if sum(rxn_1_side.values()) != sum(rxn_2_side.values()):
        logger.warning(
            '      + Number of compounds different in the template reaction and the transformation to complete'
        )
        logger.warning('         |- '+rxn_name_1+': ' + str(rxn_1_side))
        logger.warning('         |- '+rxn_name_2+': ' + str(rxn_2_side))
        # exit()
        return False
    return True

# def __round_stoichio(
#     value: float,
#     cmpd_id: str,
#     logger: Logger = getLogger(__file__)
# ) -> int:
#     r_value = round(value)
#     if r_value != value:
#         logger.warning(f'      + Stoichometric coefficient for {cmpd_id} ({value}) has been rounded.')
#     return r_value


def build_final_transfo(
    trans_input: Dict,
    missing_compounds: Dict,
    logger: Logger = getLogger(__name__)
) -> Dict:

    logger.debug(f'trans_input: {trans_input}')
    logger.debug(f'added_cmpds: {missing_compounds}')
    compl_transfo = {}

    # Add compounds to add to input transformation
    for side in Reaction.get_SIDES():

        compl_transfo[side] = deepcopy(trans_input[side])

        # All infos (stoichio, cid, smiles, InChI...) for compounds with known structures
        for cmpd_id, cmpd_sto in missing_compounds[side].items():
            if cmpd_id not in compl_transfo[side]:
                compl_transfo[side][cmpd_id] = 0
            compl_transfo[side][cmpd_id] += cmpd_sto

    logger.debug('COMPLETED TRANSFORMATION:'+str(dumps(compl_transfo, indent=4)))

    return compl_transfo


# def build_final_transfo_legacy(
#     trans_input: Dict,
#     missing_compounds: Dict,
#     logger: Logger = getLogger(__name__)
# ) -> Dict:

#     logger.debug(f'trans_input: {trans_input}')
#     logger.debug(f'added_cmpds: {missing_compounds}')
#     compl_transfo = {}

#     # Add compounds to add to input transformation
#     for side in Reaction.get_SIDES():

#         compl_transfo[side] = deepcopy(trans_input[side])

#         # All infos (stoichio, cid, smiles, InChI...) for compounds with known structures
#         for cmpd_id, cmpd_infos in missing_compounds[side].items():
#             _cmpd_id = cmpd_infos[trans_input['format']]
#             # r_stoichio = __round_stoichio(cmpd_infos['stoichio'], cmpd_id, logger)
#             if _cmpd_id not in compl_transfo[side]:
#                 compl_transfo[side][_cmpd_id] = 0
#             compl_transfo[side][_cmpd_id] += cmpd_infos['stoichio']

#         # Only cid and stoichio for compounds with no structure
#         for cmpd_id, cmpd_infos in missing_compounds[side+'_nostruct'].items():
#             _cmpd_id = cmpd_infos['cid']
#             # r_stoichio = __round_stoichio(cmpd_infos['stoichio'], cmpd_id, logger)
#             if _cmpd_id not in compl_transfo[side]:
#                 compl_transfo[side][_cmpd_id] = 0
#             compl_transfo[side][_cmpd_id] += cmpd_infos['stoichio']

#     logger.debug('COMPLETED TRANSFORMATION:'+str(dumps(compl_transfo, indent=4)))

#     return compl_transfo


def build_trans_input(
    transfo: str,
    logger: Logger = getLogger(__file__)
) -> Dict:
    """
    Build transformation to complete.

    Parameters
    ----------
    trans_smi: str
        Transformation in SMILES format or with CID.
        Stoichiometric coefficients must be separated by spaces.
    logger : Logger
        The logger object.

    Returns
    -------
    transfo: Dict
        Dictionary of the transformation.
    """
    logger.debug(f'transfo: {transfo}')

    trans_input = {
        'left': {},
        'right': {},
        'format': '',
        'sep_side': '',
        'sep_cmpd': ''
    }
    # Detect input format
    if '>>' in transfo:  # SMILES
        trans_input['format'] = 'smiles'
        trans_input['sep_side'] = '>>'
        trans_input['sep_cmpd'] = '.'
    elif '=' in transfo:  # CMPD IDs
        trans_input['format'] = 'cid'
        trans_input['sep_side'] = '='
        trans_input['sep_cmpd'] = '+'
    trans = {}
    trans['left'], trans['right'] = transfo.split(trans_input['sep_side'])
    for side in Reaction.get_SIDES():
        for cmpd in trans[side].split(trans_input['sep_cmpd']):
            # Separate compounds, remove leading and trailing spaces
            _list = cmpd.strip().split(' ')
            # Detect stoichio coeff
            if len(_list) > 1:
                _coeff = float(_list[0])
                _cmpd = _list[1]
            else:
                _coeff = 1.0
                _cmpd = _list[0]
            if not _cmpd in trans_input[side]:
                trans_input[side][_cmpd] = 0
            trans_input[side][_cmpd] += _coeff
    logger.debug('INPUT TRANSFORMATION: '+str(dumps(trans_input, indent=4)))
    return trans_input


def detect_missing_compounds(
    tmpl_rxn: Dict,
    rxn_rule: Dict,
    cid_strc: Dict,
    cmpds_to_ignore: List[str] = [],
    logger: Logger = getLogger(__file__)
) -> Tuple[Dict, List]:
    """
    Find compounds to be added from template reaction and reaction rule.

    Parameters
    ----------
    tmpl_rxn: Dict
        Orignial reaction.
    rxn_rule: Dict
        Reaction rule.
    cid_strc: Dict
        Compound structures.
    cmpds_to_ignore: List[str]
        List of compounds to ignore.
    logger : Logger
        The logger object.

    Returns
    -------
    added_compounds: Dict
        Compounds to add in various formats
    compounds_nostruct: List
        Compounds with no structure
    """

    logger.debug(f'tmpl_rxn: {tmpl_rxn}')
    logger.debug(f'rxn_rule: {rxn_rule}')

    added_compounds = {
        'left': {},
        'right': {},
        'left_nostruct': {},
        'right_nostruct': {}
    }

    for side in Reaction.get_SIDES():
        # Get the difference between the template reaction and the reaction rule,
        # difference in compounds and in stoichio coeff
        diff_cmpds = dict(
            (Counter(tmpl_rxn[side]) - Counter(rxn_rule[side]))
        )
        # Fill the dictionary with all informations about the compounds to add
        for cmp_id, cmp_sto in diff_cmpds.items():
            if cmp_id in cmpds_to_ignore:
                logger.warning(f'      + Ignoring compound {cmp_id} ({cmp_sto}) on {side} side of the transformation')
                continue
            # Handle compounds with no structure
            if cmp_id in cid_strc and cid_strc[cmp_id]['smiles'] not in [None, '']:
                added_compounds[side][cmp_id] = cmp_sto
            else:
                logger.warning(f'      Compound {cmp_id} with no structure not added to the transformation. Please check the template reaction and the reaction rule, and update the cache if needed.')
                # Fill only 'stoichio' information
                # added_compounds[side+'_nostruct'][cmp_id] = cmp_sto
                ## added_compounds[side+'_nostruct'][cmp_id] = {}
                ## added_compounds[side+'_nostruct'][cmp_id]['stoichio'] = cmp_sto
                ## added_compounds[side+'_nostruct'][cmp_id]['cid'] = cmp_id

    return added_compounds


# def detect_missing_compounds(
#     rxn: Dict,
#     rxn_ref: Dict,
#     cid_strc: Dict,
#     cmpds_to_ignore: List[str] = [],
#     logger: Logger = getLogger(__file__)
# ) -> Tuple[Dict, List]:
#     """
#     Find compounds to be added to a reaction from a reference reaction.

#     Parameters
#     ----------
#     rxn: Dict
#         Reaction to detect missing compounds.
#     rxn_ref: Dict
#         Reference rule.
#     cid_strc: Dict
#         Compound structures.
#     cmpds_to_ignore: List[str]
#         List of compounds to ignore.
#     logger : Logger
#         The logger object.

#     Returns
#     -------
#     added_compounds: Dict
#         Compounds to add in various formats
#     compounds_nostruct: List
#         Compounds with no structure
#     """

#     logger.debug(f'rxn: {rxn}')
#     logger.debug(f'rxn_ref: {rxn_ref}')

#     added_compounds = {
#         'left': {},
#         'right': {},
#         'left_nostruct': {},
#         'right_nostruct': {}
#     }

#     for side in Reaction.get_SIDES():
#         # Get the difference with the reference reaction,
#         # difference in compounds and in stoichio coeff
#         diff_cmpds = dict(
#             (Counter(rxn_ref[side]) - Counter(rxn[side]))
#         )
#         logger.debug(f'Reference reaction {side} side: {rxn_ref[side]}')
#         logger.debug(f'Reaction {side} side: {rxn[side]}')
#         logger.debug(f'DIFFERENCE on {side} side: {diff_cmpds}')
#         # print()
#         # print(side)
#         # print("="*len(side))
#         # print("TMPL:", tmpl_rxn[side])
#         # print("RULE:", rxn_rule[side])
#         # print("DIFF:", diff_cmpds)
#         # print()
#         # Fill the dictionary with all informations about the compounds to add
#         logger.debug(f'Detecting missing compounds on {side} side: {diff_cmpds}')
#         for cmp_id, cmp_sto in diff_cmpds.items():
#             if cmp_id in cmpds_to_ignore:
#                 logger.warning(f'      + Ignoring compound {cmp_id} ({cmp_sto}) on {side} side of the transformation')
#                 continue
#             # Handle compounds with no structure
#             if cmp_id in cid_strc and cid_strc[cmp_id]['smiles'] not in [None, '']:
#                 logger.debug('      Adding compound with structure: '+cmp_id)
#                 added_compounds[side][cmp_id] = {'stoichio': cmp_sto}
#                 for key, val in cid_strc[cmp_id].items():
#                     # add val from key
#                     if key == 'mnxm':
#                         key = 'cid'
#                     added_compounds[side][cmp_id][key] = val
#             else:
#                 logger.warning('      Compounds with no structure: '+cmp_id)
#                 # Fill only 'stoichio' information
#                 added_compounds[side+'_nostruct'][cmp_id] = {'stoichio': cmp_sto, 'cid': cmp_id}

#     logger.debug('ADDED COMPOUNDS: '+str(dumps(added_compounds, indent=4)))

#     return added_compounds
