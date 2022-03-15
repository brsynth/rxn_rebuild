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


SIDES = ['left', 'right']


def rebuild_rxn(
    rxn_rule_id: str,
    transfo: str,
    direction: str = 'reverse',
    tmpl_rxn_id: str = None,
    cache: 'rrCache' = None,
    logger: Logger = getLogger(__name__)
) -> str:

    logger.debug(f'rxn_rule_id: {rxn_rule_id}')
    logger.debug(f'transfo: {transfo}')
    logger.debug(f'tmpl_rxn_id: {tmpl_rxn_id}')
    logger.debug(f'direction: {direction}')

    ## INPUT TRANSFORMATION
    trans_input = build_trans_input(transfo, logger)  # remove whitespaces
    # # If input transfo is provided in forward direction,
    # # swap left and right sides
    # if direction == 'forward' or direction == 'fwd':
    #     trans_input['left'], trans_input['right'] = trans_input['right'], trans_input['left']

    ## LOAD CACHE
    if cache is None:
        cache = rrCache(
            attrs=['rr_reactions', 'template_reactions', 'cid_strc']
            # logger=logger
        )

    ## COMPLETE TRANSFORMATION
    completed_transfos = {}
    if not tmpl_rxn_id is None:
        completed_transfos[tmpl_rxn_id] = complete_transfo(
            trans_input=trans_input,
            direction=direction,
            rxn_rule=cache.get('rr_reactions')[rxn_rule_id][tmpl_rxn_id],
            tmpl_rxn=cache.get('template_reactions')[tmpl_rxn_id],
            tmpl_rxn_id=tmpl_rxn_id,
            compounds=cache.get('cid_strc'),
            logger=logger
        )
    else:  # One completed transformation per template reaction
        for tpl_rxn_id, rxn_rule in cache.get('rr_reactions')[rxn_rule_id].items():
            completed_transfos[tpl_rxn_id] = complete_transfo(
                trans_input=trans_input,
                direction=direction,
                rxn_rule=rxn_rule,
                tmpl_rxn=cache.get('template_reactions')[tpl_rxn_id],
                tmpl_rxn_id=tpl_rxn_id,
                compounds=cache.get('cid_strc'),
                logger=logger
            )

    return completed_transfos


def complete_transfo(
    trans_input: Dict,
    rxn_rule: Dict,
    tmpl_rxn: Dict,
    tmpl_rxn_id: str,
    compounds: Dict,
    direction: str = 'reverse',
    logger: Logger = getLogger(__name__)
) -> Dict:

    logger.debug('REACTION RULE:'+str(dumps(rxn_rule, indent=4)))

    # If input transfo is provided in forward direction,
    # swap left and right sides of reaction rule, and
    # change 'rel_direction'
    if direction == 'forward' or direction == 'fwd':
        _rxn_rule = {
            'left': rxn_rule['right'],
            'right': rxn_rule['left'],
            'rel_direction': -rxn_rule['rel_direction']
        }
    else:  # reverse
        _rxn_rule = deepcopy(rxn_rule)

    ## CHECK 1/2
    # Check if the number of structures in the right part of SMILES of transformation to complete
    # is equal to the number of products of the template reaction used to build the reaction rule.
    # Just in right part since rules are always mono-component
    check_compounds_number(
        'INPUT TRANSFORMATION',
        trans_input,
        'REACTION RULE',
        _rxn_rule,
        'right',
        logger=logger
    )

    ## TEMPLATE REACTION
    if _rxn_rule['rel_direction'] == -1:
        _tmpl_rxn = {
            'left': tmpl_rxn['right'],
            'right': tmpl_rxn['left'],
        }
    else:
        _tmpl_rxn = deepcopy(tmpl_rxn)

    ## ADD MISSING COMPOUNDS TO THE FINAL TRANSFORMATION
    missing_compounds = detect_missing_compounds(
        _tmpl_rxn,
        _rxn_rule,
        compounds,
        logger=logger
    )

    ## BUILD FINAL TRANSFORMATION
    compl_transfo = build_final_transfo(
        trans_input=trans_input,
        missing_compounds=missing_compounds,
        logger=logger
    )

    ## CHECK 2/2
    # Check if the number of compounds in both right and left sides of SMILES of the completed transformation
    # is equal to the ones of the template reaction
    for side in SIDES:
        check_compounds_number(
            'COMPLETED TRANSFORMATION',
            compl_transfo,
            'TEMPLATE REACTION ({tmpl_rxn_id})'.format(tmpl_rxn_id=tmpl_rxn_id),
            _tmpl_rxn,
            side,
            logger=logger
        )

    return {
        'full_transfo': compl_transfo,
        'added_cmpds': missing_compounds,
        'sep_side': trans_input['sep_side'],
        'sep_cmpd': trans_input['sep_cmpd']
    }


def check_compounds_number(
    rxn_name_1: str,
    rxn_1: Dict,
    rxn_name_2: str,
    rxn_2: Dict,
    side: str,
    logger: Logger = getLogger(__name__)
) -> None:
    logger.debug(f'rxn_1: {rxn_1}')
    logger.debug(f'rxn_2: {rxn_2}')
    # Check if the number of structures in the part of SMILES of rxn_1
    # is equal to the number of products of rxn_2.
    if sum(rxn_1[side].values()) != sum(rxn_2[side].values()):
        logger.warning(
            '      + Number of compounds different in the template reaction and the transformation to complete'
        )
        logger.warning('         |- '+rxn_name_1+' ['+side+']: ' + str(rxn_1[side]))
        logger.warning('         |- '+rxn_name_2+' ['+side+']: ' + str(rxn_2[side]))


def __round_stoichio(
    value: float,
    cmpd_id: str,
    logger: Logger = getLogger(__file__)
) -> int:
    r_value = round(value)
    if r_value != value:
        logger.warning(f'      + Stoichometric coefficient for {cmpd_id} ({value}) has been rounded.')
    return r_value


def build_final_transfo(
    trans_input: Dict,
    missing_compounds: Dict,
    logger: Logger = getLogger(__name__)
) -> Dict:
    logger.debug(f'trans_input: {trans_input}')
    logger.debug(f'added_cmpds: {missing_compounds}')
    compl_transfo = {}
    # Add compounds to add to input transformation
    for side in SIDES:
        compl_transfo[side] = deepcopy(trans_input[side])
        # All infos (stoichio, cid, smiles, InChI...) for compounds with known structures a
        for cmpd_id, cmpd_infos in missing_compounds[side].items():
            _cmpd_id = cmpd_infos[trans_input['format']]
            # r_stoichio = __round_stoichio(cmpd_infos['stoichio'], cmpd_id, logger)
            if not _cmpd_id in compl_transfo[side]:
                compl_transfo[side][_cmpd_id] = 0
            compl_transfo[side][_cmpd_id] += cmpd_infos['stoichio']
        # Only cid and stoichio for compounds with no structure
        for cmpd_id, cmpd_infos in missing_compounds[side+'_nostruct'].items():
            _cmpd_id = cmpd_infos['cid']
            # r_stoichio = __round_stoichio(cmpd_infos['stoichio'], cmpd_id, logger)
            if not _cmpd_id in compl_transfo[side]:
                compl_transfo[side][_cmpd_id] = 0
            compl_transfo[side][_cmpd_id] += cmpd_infos['stoichio']
    logger.debug('COMPLETED TRANSFORMATION:'+str(dumps(compl_transfo, indent=4)))
    return compl_transfo


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
    for side in SIDES:
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

    for side in SIDES:
        # Get the difference between the template reaction and the reaction rule,
        # difference in compounds and in stoichio coeff
        diff_cmpds = dict(
            (Counter(tmpl_rxn[side]) - Counter(rxn_rule[side]))
        )
        # print()
        # print(side)
        # print("="*len(side))
        # print("TMPL:", tmpl_rxn[side])
        # print("RULE:", rxn_rule[side])
        # print("DIFF:", diff_cmpds)
        # print()
        # Fill the dictionary with all informations about the compounds to add
        for cmp_id, cmp_sto in diff_cmpds.items():
            # Handle compounds with no structure
            if cmp_id in cid_strc:
                added_compounds[side][cmp_id] = {}
                added_compounds[side][cmp_id]['stoichio'] = cmp_sto
                for key, val in cid_strc[cmp_id].items():
                    # add val from key
                    if key == 'mnxm':
                        key = 'cid'
                    added_compounds[side][cmp_id][key] = val
            else:
                logger.warning('      Compounds with no structure: '+cmp_id)
                # Fill only 'stoichio' information
                added_compounds[side+'_nostruct'][cmp_id] = {}
                added_compounds[side+'_nostruct'][cmp_id]['stoichio'] = cmp_sto
                added_compounds[side+'_nostruct'][cmp_id]['cid'] = cmp_id

    return added_compounds
