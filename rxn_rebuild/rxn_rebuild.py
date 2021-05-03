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

    ## REACTION RULE
    rxn_rules = load_rxn_rules(
        cache,
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
        tmpl_rxn = load_tmpl_rxn(cache, tmpl_rxn_id, rxn_rule['rel_direction'])
        logger.debug('TEMPLATE REACTION:'+str(dumps(tmpl_rxn, indent=4)))

        completed_transfos[tmpl_rxn_id] = {}

        ## ADD MISSING COMPOUNDS TO THE FINAL TRANSFORMATION
        added_compounds = add_compounds(
            cache,
            tmpl_rxn,
            rxn_rule
        )
        completed_transfos[tmpl_rxn_id]['added_cmpds'] = added_compounds

        ## BUILD FINAL TRANSFORMATION
        compl_transfo = {}
        # Add compound to add to input transformation
        for side in ['left', 'right']:
            compl_transfo[side] = list(trans_input[side])
            for cmpd_id, cmpd_infos in completed_transfos[tmpl_rxn_id]['added_cmpds'][side].items():
                compl_transfo[side] += [cmpd_infos[trans_input['format']]]*cmpd_infos['stoichio']
        logger.debug('COMPLETED TRANSFORMATION:'+str(dumps(compl_transfo, indent=4)))

        completed_transfos[tmpl_rxn_id]['full_transfo'] = trans_input['sep_cmpd'].join(compl_transfo['left'])+trans_input['sep_side']+trans_input['sep_cmpd'].join(compl_transfo['right'])

        # Check if the number of compounds in both right and left sides of SMILES of the completed transformation
        # is equal to the ones of the template reaction
        for side in ['left', 'right']:
            if len(compl_transfo[side]) != sum(tmpl_rxn[side].values()):
                logger.warning(
                    '      Number of compounds different in the template reaction and the completed transformation'
                )
                logger.warning('         |- COMPLETED TRANSFORMATION ['+side+']: ' + str(compl_transfo[side]))
                logger.warning('         |- TEMPLATE REACTION ['+side+']: ' + str(tmpl_rxn[side]))

    return completed_transfos


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
        trans_input['format'] = 'id'
        trans_input['sep_side'] = '='
        trans_input['sep_cmpd'] = '+'
    trans_left, trans_right = transfo.split(trans_input['sep_side'])
    for cmpd in trans_left.split(trans_input['sep_cmpd']):
        trans_input['left'] += [cmpd]
    for cmpd in trans_right.split(trans_input['sep_cmpd']):
        trans_input['right'] += [cmpd]
    return trans_input


def load_rxn_rules(
    cache: 'rrCache',
    rxn_rule_id: str,
    tmpl_rxn_id: str,
    logger: Logger=getLogger(__file__)
) -> Tuple[str, Dict]:
    """
    Seeks for the reaction rule of ID rxn_rule_id in the cache.

    Parameters
    ----------
    cache: rrCache
        Pre-computed data.
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
    cache.load(['rr_reactions'])
    rxn_rules = {}
    # Get the reaction rule built from the template reaction of ID 'tmpl_rxn_id'
    if tmpl_rxn_id is not None:
        rxn_rules[tmpl_rxn_id] = cache.get('rr_reactions')[rxn_rule_id][tmpl_rxn_id]
    else: # Get all template reactions
        for rxn_id in cache.get('rr_reactions')[rxn_rule_id].keys():
            rxn_rules[rxn_id] = cache.get('rr_reactions')[rxn_rule_id][rxn_id]
    return rxn_rules


def load_tmpl_rxn(
    cache: 'rrCache',
    rxn_id: str,
    rr_direction: int,
    logger: Logger=getLogger(__file__)
) -> Dict:
    """
    Seeks for the template reaction of ID rxn_id in the cache.

    Parameters
    ----------
    cache: rrCache
        Pre-computed data.
    rxn_id: str
        ID of the orignial reaction.
    rr_direction: int
        Direction of the reaction used to build the rule.
    logger : Logger
        The logger object.

    Returns
    -------
    tmpl_rxn: Dict
        template reaction looked for.
    """
    cache.load(['rr_full_reactions'])
    rxn = cache.get('rr_full_reactions')[rxn_id]
    if rr_direction == 1:
        rxn_left  = rxn['left']
        rxn_right = rxn['right']
    else:
        left = dict(rxn['left'])
        rxn['left']  = rxn['right']
        rxn['right'] = left
    return rxn


def add_compounds(
    cache: 'rrCache',
    tmpl_rxn: Dict,
    rxn_rule: Dict,
    logger: Logger=getLogger(__file__)
) -> Dict:
    """
    Find compounds to be added from template reaction and reaction rule.

    Parameters
    ----------
    cache: rrCache
        Pre-computed data.
    tmpl_rxn: Dict
        Orignial reaction.
    rxn_rule: Dict
        Reaction rule.
    logger : Logger
        The logger object.

    Returns
    -------
    trans_res: Dict
        Compounds to add in various formats
    """
    added_compounds = {
        'left': {},
        'right': {}
        # 'smiles': {'left': {}, 'right': {}},
        # 'id': {'left': {}, 'right': {}},
        # 'formula': {'left': {}, 'right': {}},
        # 'inchi': {'left': {}, 'right': {}},
        # 'inchikey': {'left': {}, 'right': {}},
        # 'name': {'left': {}, 'right': {}}
    }
    cache.load(['cid_strc'])

    # tmpl_rxn['right']['MNXM821'] = 3

    for side in ['left', 'right']:
        # Get the difference between the template reaction and the reaction rule,
        # difference in compounds and in stoichio coeff
        diff_cmpds = dict(
            (Counter(tmpl_rxn[side]) - Counter(rxn_rule[side]))
        )
        # Fill the dictionary with all informations about the compounds to add
        for cmp_id, cmp_sto in diff_cmpds.items():
            added_compounds[side][cmp_id] = {}
            added_compounds[side][cmp_id]['stoichio'] = cmp_sto
            for key, val in cache.get('cid_strc')[cmp_id].items():
                # add val from key
                if key == 'mnxm':
                    key = 'id'
                added_compounds[side][cmp_id][key] = val

    return added_compounds

    # ## LEFT
    # for cmp_id in detected_compounds['toadd']['left']:
    #     # get smiles from compound ID
    #     smi = cache.get('cid_strc')[cmp_id]['smiles']
    #     # add the compound x sto coeff in the template reaction
    #     trans_res['smiles']['left'].extend([smi]*tmpl_rxn['left'][cmp_id])
    # # add elements already in input transfo because sto coeff could be > 1 in template reaction
    # for cmp_id in detected_compounds['common']['left']:
    #     # get smiles from compound ID
    #     smi = cache.get('cid_strc')[cmp_id]['smiles']
    #     # add the compound x-1 (already in input transfo) sto coeff in the template reaction
    #     trans_res['smiles']['left'].extend([smi]*(tmpl_rxn['left'][cmp_id]-1))

    # ## RIGHT
    # for cmp_id in detected_compounds['toadd']['right']:
    #     # get smiles from compound ID
    #     smi = cache.get('cid_strc')[cmp_id]['smiles']
    #     # add the compound x sto coeff in the template reaction
    #     trans_res['smiles']['right'].extend([cmp_id]*tmpl_rxn['right'][cmp_id])
    # # add elements already in input transfo because sto coeff could be > 1 in template reaction
    # for cmp_id in detected_compounds['common']['right']:
    #     # get smiles from compound ID
    #     smi = cache.get('cid_strc')[cmp_id]['smiles']
    #     # add the compound x-1 (already in input transfo) sto coeff in the template reaction
    #     trans_res['smiles']['right'].extend([cmp_id]*(tmpl_rxn['right'][cmp_id]-1))

    # return trans_res['smiles']