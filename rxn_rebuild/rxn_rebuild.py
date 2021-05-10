from logging import (
    Logger,
    getLogger,
    StreamHandler
)
from typing import (
    List,
    Dict,
    Tuple
)
from rr_cache import rrCache
from brs_utils import diff
from collections import Counter
from json import dumps


SIDES = ['left', 'right']


def rebuild_rxn(
    rxn_rule_id: str,
    transfo: str,
    tmpl_rxn_id: str = None,
    cache: 'rrCache'=None,
    logger: Logger=getLogger(__name__)
) -> str:

    ## INPUT TRANSFORMATION
    trans_input = build_trans_input(transfo.replace(' ', '')) # remove whitespaces

    ## LOAD CACHE
    if cache is None:
        cache = rrCache(
            db='file',
            attrs=['rr_reactions', 'template_reactions','cid_strc']
            # logger=logger
        )

    ## COMPLETE TRANSFORMATION
    completed_transfos = {}
    if tmpl_rxn_id is not None:
        completed_transfos[tmpl_rxn_id] = complete_transfo(
            trans_input,
            cache.get('rr_reactions')[rxn_rule_id][tmpl_rxn_id],
            tmpl_rxn_id,
            cache,
            logger=logger
        )
    else: # One completed transformation per template reaction
        for tpl_rxn_id, rxn_rule in cache.get('rr_reactions')[rxn_rule_id].items():
            completed_transfos[tpl_rxn_id] = complete_transfo(
                trans_input,
                rxn_rule,
                tpl_rxn_id,
                cache,
                logger=logger
            )

    return completed_transfos


def complete_transfo(
    trans_input: Dict,
    rxn_rule: Dict,
    tmpl_rxn_id: str,
    cache: 'rrCache',
    logger: Logger=getLogger(__name__)
) -> Dict:

    completed_transfo = {}

    logger.debug('REACTION RULE:'+str(dumps(rxn_rule, indent=4)))

    ## CHECK 1/2
    # Check if the number of structures in the right part of SMILES of transformation to complete
    # is equal to the number of products of the template reaction used to build the reaction rule.
    # Just in right part since rules are always mono-component
    check_compounds_number(
        'INPUT TRANSFORMATION',
        trans_input,
        'REACTION RULE',
        rxn_rule,
        'right',
        logger=logger
    )

    ## TEMPLATE REACTION
    tmpl_rxn = load_tmpl_rxn(
        cache.get('template_reactions'),
        tmpl_rxn_id,
        rxn_rule['rel_direction'],
        logger=logger
    )

    ## ADD MISSING COMPOUNDS TO THE FINAL TRANSFORMATION
    missing_compounds = detect_missing_compounds(
        tmpl_rxn,
        rxn_rule,
        cache.get('cid_strc'),
        logger=logger
    )

    ## BUILD FINAL TRANSFORMATION
    compl_transfo = build_final_transfo(
        trans_input,
        missing_compounds,
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
            tmpl_rxn,
            side,
            logger=logger
        )

    # Build full transformation only if there is no compound without structure
    # or if the input format of the transforamtion to complete is 'cid'
    if missing_compounds['left_nostruct'] == missing_compounds['right_nostruct'] == {} \
    or trans_input['format'] == 'cid':
        completed_transfo['full_transfo'] = \
              trans_input['sep_cmpd'].join(compl_transfo['left']) \
            + trans_input['sep_side'] \
            + trans_input['sep_cmpd'].join(compl_transfo['right'])

    completed_transfo['added_cmpds'] = missing_compounds

    return completed_transfo

def check_compounds_number(
    rxn_name_1: str,
    rxn_1: Dict,
    rxn_name_2: str,
    rxn_2: Dict,
    side: str,
    logger: Logger=getLogger(__name__)
) -> None:
    # Check if the number of structures in the part of SMILES of rxn_1
    # is equal to the number of products of rxn_2.
    if len(rxn_1[side]) != sum(rxn_2[side].values()):
        logger.warning(
            '      Number of compounds different in the template reaction and the transformation to complete'
        )
        logger.warning('         |- '+rxn_name_1+' ['+side+']: ' + str(rxn_1[side]))
        logger.warning('         |- '+rxn_name_2+' ['+side+']: ' + str(rxn_2[side]))


def build_final_transfo(
    trans_input: Dict,
    added_cmpds: Dict,
    logger: Logger=getLogger(__name__)
) -> Dict:
    compl_transfo = {}
    # Add compound to add to input transformation
    for side in SIDES:
        compl_transfo[side] = list(trans_input[side])
        # All infos (stoichio, cid, smiles, InChI...) for compounds with known structures a
        for cmpd_id, cmpd_infos in added_cmpds[side].items():
            compl_transfo[side] += [cmpd_infos[trans_input['format']]]*cmpd_infos['stoichio']
        # Only cid and stoichio for compounds with no structure
        for cmpd_id, cmpd_infos in added_cmpds[side+'_nostruct'].items():
            compl_transfo[side] += [cmpd_infos['cid']]*cmpd_infos['stoichio']
    logger.debug('COMPLETED TRANSFORMATION:'+str(dumps(compl_transfo, indent=4)))
    return compl_transfo


def build_trans_input(
    transfo: str,
    logger: Logger=getLogger(__file__)
) -> Dict:
    """
    Build transformation to complete.

    Parameters
    ----------
    trans_smi: str
        Transformation in SMILES format.
    logger : Logger
        The logger object.

    Returns
    -------
    transfo: Dict
        Dictionary of the transformation.
    """
    trans_input = {
        'left': [],
        'right': [],
        'format': '',
        'sep_side': '',
        'sep_cmpd': ''
    }
    # Detect input format
    if '>>' in transfo: # SMILES
        trans_input['format'] = 'smiles'
        trans_input['sep_side'] = '>>'
        trans_input['sep_cmpd'] = '.'
    elif '=' in transfo: # CMPD IDs
        trans_input['format'] = 'cid'
        trans_input['sep_side'] = '='
        trans_input['sep_cmpd'] = '+'
    trans = {}
    trans['left'], trans['right'] = transfo.split(trans_input['sep_side'])
    for side in SIDES:
        for cmpd in trans[side].split(trans_input['sep_cmpd']):
            trans_input[side] += [cmpd]
    logger.debug('INPUT TRANSFORMATION: '+str(dumps(trans_input, indent=4)))
    return trans_input


def load_tmpl_rxn(
    tmplt_rxns: Dict,
    rxn_id: str,
    rr_direction: int,
    logger: Logger=getLogger(__file__)
) -> Dict:
    """
    Seeks for the template reaction of ID rxn_id in the cache.

    Parameters
    ----------
    tmplt_rxns: Dict
        Template (original) reactions.
    rxn_id: str
        ID of the template reaction.
    rr_direction: int
        Direction of the reaction used to build the rule.
    logger : Logger
        The logger object.

    Returns
    -------
    tmpl_rxn: Dict
        template reaction looked for.
    """
    rxn = tmplt_rxns[rxn_id]
    if rr_direction == 1:
        rxn_left  = rxn['left']
        rxn_right = rxn['right']
    else:
        left = dict(rxn['left'])
        rxn['left']  = rxn['right']
        rxn['right'] = left
    logger.debug('TEMPLATE REACTION:'+str(dumps(rxn, indent=4)))
    return rxn


def detect_missing_compounds(
    tmpl_rxn: Dict,
    rxn_rule: Dict,
    cid_strc: Dict,
    logger: Logger=getLogger(__file__)
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
    added_compounds = {
        'left': {},
        'right': {},
        'left_nostruct': {},
        'right_nostruct': {}
    }
    compounds_nostruct = []

    for side in SIDES:
        # Get the difference between the template reaction and the reaction rule,
        # difference in compounds and in stoichio coeff
        diff_cmpds = dict(
            (Counter(tmpl_rxn[side]) - Counter(rxn_rule[side]))
        )
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
