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
    trans_smi: str,
    tmpl_rxn_id: str = None,
    logger: Logger=getLogger(__name__)
) -> str:

    ## INPUT TRANSFORMATION
    trans_input = build_trans_input(trans_smi)
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

        ## ADD MISSING COMPOUNDS INTO FINAL TRANSFORMATION
        trans_res = complete_input_trans(
            cache,
            trans_input,
            tmpl_rxn,
            rxn_rule
        )
        logger.debug('COMPLETED TRANSFORMATION:'+str(dumps(trans_res, indent=4)))

        # Check if the number of compounds in both right and left parts of SMILES of the completed transformation
        # is equal to the ones of the template reaction
        for side in ['left', 'right']:
            if len(trans_res[side]) != sum(tmpl_rxn[side].values()):
                logger.warning(
                    '      Number of compounds different in the template reaction and the completed transformation'
                )
                logger.warning('         |- COMPLETED TRANSFORMATION ['+side+']: ' + str(trans_res[side]))
                logger.warning('         |- TEMPLATE REACTION ['+side+']: ' + str(tmpl_rxn[side]))

        completed_transfos[tmpl_rxn_id] = '.'.join(trans_res['left'])+'>>'+'.'.join(trans_res['right'])

    return completed_transfos


def build_trans_input(
    trans_smi: str,
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
        'right': []
    }
    trans_left, trans_right = trans_smi.split('>>')
    for smi in trans_left.split('.'):
        trans_input['left'] += [smi]
    for smi in trans_right.split('.'):
        trans_input['right'] += [smi]
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
    if tmpl_rxn_id is not None:
        rxn_rules[tmpl_rxn_id] = cache.get('rr_reactions')[rxn_rule_id][tmpl_rxn_id]
    else: # Take all template reactions
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


def complete_input_trans(
    cache: 'rrCache',
    trans_input: Dict,
    tmpl_rxn: Dict,
    rxn_rule: Dict,
    logger: Logger=getLogger(__file__)
) -> Dict:
    """
    Complete input transformation from template reaction and reaction rule.

    Parameters
    ----------
    cache: rrCache
        Pre-computed data.
    trans_input: Dict
        Transformation to complete
    tmpl_rxn: Dict
        Orignial reaction.
    rxn_rule: Dict
        Reaction rule.
    logger : Logger
        The logger object.

    Returns
    -------
    trans_res: Dict
        Completed transformation
    """
    trans_res = {
        'smiles': {},
        'ids': {}
    }
    trans_res['smiles'] = dict(trans_input)
    cache.load(['cid_strc'])

    for side in ['left', 'right']:
        # Get the difference between the template reaction and the reaction rule,
        # difference in compounds and in stoichio coeff
        diff_cmpds = list(
            (Counter(tmpl_rxn[side]) - Counter(rxn_rule[side])).elements()
        )
        for cmp_id in diff_cmpds:
            # get smiles from compound ID
            smi = cache.get('cid_strc')[cmp_id]['smiles']
            # add the compound x sto coeff in the template reaction
            trans_res['smiles'][side] += [smi]

    return trans_res['smiles']

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