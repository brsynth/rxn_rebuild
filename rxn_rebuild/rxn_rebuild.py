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


def rebuild_rxn(
    cache: 'rrCache',
    rxn_rule_id: str,
    transfo: str,
    tmpl_rxn_id: str = None,
    logger: Logger=getLogger(__name__)
) -> str:

    ## INPUT TRANSFORMATION
    trans_input = build_trans_input(transfo.replace(' ', '')) # remove whitespaces
    logger.debug('INPUT TRANSFORMATION: '+str(dumps(trans_input, indent=4)))

    ## LOAD CACHE
    load_cache(cache)

    ## REACTION RULES
    rxn_rules = load_rxn_rules(
        cache.get('rr_reactions'),
        rxn_rule_id,
        tmpl_rxn_id
    )

    completed_transfos = {}

    for tmpl_rxn_id, rxn_rule in rxn_rules.items():
        logger.debug('REACTION RULE:'+str(dumps(rxn_rule, indent=4)))

        # Check if the number of structures in the right part of SMILES of transformation to complete
        # is equal to the number of products of the template reaction used to build the reaction rule.
        # Just in right part since rules are always mono-component
        if len(trans_input['right']) != sum(rxn_rule['right'].values()):
            logger.warning(
                '      Number of compounds different in the template reaction and the transformation to complete'
            )
            logger.warning('         |- INPUT TRANSFORMATION [right]: ' + str(trans_input['right']))
            logger.warning('         |- REACTION RULE [right]: ' + str(rxn_rule['right']))

        ## TEMPLATE REACTION
        tmpl_rxn = load_tmpl_rxn(
            cache.get('rr_full_reactions'),
            tmpl_rxn_id,
            rxn_rule['rel_direction']
        )
        logger.debug('TEMPLATE REACTION:'+str(dumps(tmpl_rxn, indent=4)))

        completed_transfos[tmpl_rxn_id] = {}

        ## ADD MISSING COMPOUNDS TO THE FINAL TRANSFORMATION
        added_compounds = add_compounds(
            cache.get('cid_strc'),
            tmpl_rxn,
            rxn_rule
        )

        completed_transfos[tmpl_rxn_id]['added_cmpds'] = added_compounds

        ## BUILD FINAL TRANSFORMATION
        compl_transfo = complete_transfo(
            trans_input,
            completed_transfos[tmpl_rxn_id]
        )
        logger.debug('COMPLETED TRANSFORMATION:'+str(dumps(compl_transfo, indent=4)))

        # Check if the number of compounds in both right and left sides of SMILES of the completed transformation
        # is equal to the ones of the template reaction
        for side in ['left', 'right']:
            if len(compl_transfo[side]) != sum(tmpl_rxn[side].values()):
                logger.warning(
                    '      Number of compounds different in the template reaction and the completed transformation'
                )
                logger.warning('         |- COMPLETED TRANSFORMATION ['+side+']: ' + str(compl_transfo[side]))
                logger.warning('         |- TEMPLATE REACTION ['+side+']: ' + str(tmpl_rxn[side]))

        if added_compounds['left_nostruct'] == {} and added_compounds['right_nostruct'] == {} or trans_input['format'] != 'smiles':
            completed_transfos[tmpl_rxn_id]['full_transfo'] = trans_input['sep_cmpd'].join(compl_transfo['left'])+trans_input['sep_side']+trans_input['sep_cmpd'].join(compl_transfo['right'])

    return completed_transfos


def complete_transfo(
    trans_input: Dict,
    rxn: Dict,
    logger: Logger=getLogger(__name__)
) -> Dict:
    compl_transfo = {}
    # Add compound to add to input transformation
    for side in ['left', 'right']:
        compl_transfo[side] = list(trans_input[side])
        # All infos (stoichio, cid, smiles, InChI...) for compounds with known structures a
        for cmpd_id, cmpd_infos in rxn['added_cmpds'][side].items():
            compl_transfo[side] += [cmpd_infos[trans_input['format']]]*cmpd_infos['stoichio']
        # Only cid and stoichio for compounds with no structure
        for cmpd_id, cmpd_infos in rxn['added_cmpds'][side+'_nostruct'].items():
            compl_transfo[side] += [cmpd_infos['cid']]*cmpd_infos['stoichio']
    return compl_transfo


def load_cache(
    cache: 'rrCache',
    logger: Logger=getLogger(__file__)
) -> Dict:
    cache.load(['rr_reactions'])
    cache.load(['rr_full_reactions'])
    cache.load(['cid_strc'])


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
    trans_left, trans_right = transfo.split(trans_input['sep_side'])
    for cmpd in trans_left.split(trans_input['sep_cmpd']):
        trans_input['left'] += [cmpd]
    for cmpd in trans_right.split(trans_input['sep_cmpd']):
        trans_input['right'] += [cmpd]
    return trans_input


def load_rxn_rules(
    rxn_rules_all: Dict,
    rxn_rule_id: str,
    tmpl_rxn_id: str,
    logger: Logger=getLogger(__file__)
) -> Tuple[str, Dict]:
    """
    Seeks for the reaction rule of ID rxn_rule_id in the cache.

    Parameters
    ----------
    rxn_rules_all: Dict
        Reaction rules.
    rxn_rule_id: str
        ID of reaction rule.
    tmpl_rxn_id: str
        ID of template reaction.
    logger : Logger
        The logger object.

    Returns
    -------
    rxn_id: str
        ID of the template reaction.
    rxn_rule: Dict
        Reaction Rule looked for.
    """
    rxn_rules = {}
    # Get the reaction rule built from the template reaction of ID 'tmpl_rxn_id'
    if tmpl_rxn_id is not None:
        rxn_rules[tmpl_rxn_id] = rxn_rules_all[rxn_rule_id][tmpl_rxn_id]
    else: # Get reaction rules built from all template reactions
        for rxn_id, rxn_rule in rxn_rules_all[rxn_rule_id].items():
            rxn_rules[rxn_id] = rxn_rule
    return rxn_rules


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
    return rxn


def add_compounds(
    cid_strc: Dict,
    tmpl_rxn: Dict,
    rxn_rule: Dict,
    logger: Logger=getLogger(__file__)
) -> Tuple[Dict, List]:
    """
    Find compounds to be added from template reaction and reaction rule.

    Parameters
    ----------
    cid_strc: Dict
        Compound structures.
    tmpl_rxn: Dict
        Orignial reaction.
    rxn_rule: Dict
        Reaction rule.
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

    for side in ['left', 'right']:
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
