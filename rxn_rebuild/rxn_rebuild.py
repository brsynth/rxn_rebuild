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
from copy import deepcopy
from json import dumps
from rr_cache import rrCache


def rebuild_rxn(
    cache: 'rrCache',
    rxn_rule_id: str,
    trans_smi: str,
    ori_rxn_id: str = None,
    logger: Logger=getLogger(__name__)
) -> str:

    ## INPUT TRANSFORMATION
    trans_input = build_trans_input(trans_smi)
    logger.debug('INPUT TRANSFORMATION: '+str(dumps(trans_input, indent=4)))

    ## REACTION RULE
    tmpl_rxn_id, rxn_rule = load_rxn_rule(cache, rxn_rule_id)
    logger.debug('REACTION RULE:'+str(dumps(rxn_rule, indent=4)))

    # Check if the number of structures in the right part of SMILES of transformation to complete
    # is equal to the number of products of the template reaction used to build the reaction rule.
    # Just in right part since rules are always mono-component
    if len(trans_input['right']) != sum(rxn_rule['right'].values()):
        logger.warning(
            '      Number of compounds different in the template reaction and the input transformation'
        )
        logger.warning('         |- INPUT TRANSFORMATION [right]: ' + str(trans_input['right']))
        logger.warning('         |- TEMPLATE REACTION [right]: ' + str(rxn_rule['right']))

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
                '      Number of compounds different in the template reaction and the input transformation'
            )
            logger.warning('         |- COMPLETED TRANSFORMATION ['+side+']: ' + str(trans_res[side]))
            logger.warning('         |- TEMPLATE REACTION ['+side+']: ' + str(tmpl_rxn[side]))

    return '.'.join(trans_res['left'])+'>>'+'.'.join(trans_res['right'])


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


def load_rxn_rule(
    cache: 'rrCache',
    rxn_rule_id: str,
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
    # Take the first reaction from the list of template reactions
    rxn_id = next(iter(cache.get('rr_reactions')[rxn_rule_id]))
    return rxn_id, cache.get('rr_reactions')[rxn_rule_id][rxn_id]


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
        left = deepcopy(rxn['left'])
        rxn['left']  = rxn['right']
        rxn['right'] = left
    return rxn


def detect_compounds(
    tmpl_rxn: Dict,
    rxn_rule: Dict,
    logger: Logger=getLogger(__file__)
) -> Dict:
    """
    Detect compounds that are different and the same between the template reaction and the reaction rule.

    Parameters
    ----------
    tmpl_rxn: Dict
        Orignial reaction.
    rxn_rule: Dict
        Reaction rule.
    logger : Logger
        The logger object.

    Returns
    -------
    det_cmp: Dict
        Compounds:
            - different between template reaction and the reaction rule
            - same between template reaction and the reaction rul
    """
    # building the sets
    rxn_left_set = set(tmpl_rxn['left'].keys())
    rxn_right_set = set(tmpl_rxn['right'].keys())
    rr_left_set = set(rxn_rule['left'].keys())
    rr_right_set = set(rxn_rule['right'].keys())
    return {
        # remove identical elements
        'toadd': {
            'left': rxn_left_set - rr_left_set,
            'right': rxn_right_set - rr_right_set,
        },
        # keep identical elements to check stoichio coeff
        'common': {
            'left': rxn_left_set.intersection(rr_left_set),
            'right': rxn_right_set.intersection(rr_right_set),
        }
    }


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
    trans_res = deepcopy(trans_input)
    cache.load(['cid_strc'])

    ## DETECT IDENTICAL COMPOUNDS IN BOTH REACTION RULE AND TEMPLATE REACTION
    detected_compounds = detect_compounds(tmpl_rxn, rxn_rule)

    # LEFT
    for cmp_id in detected_compounds['toadd']['left']:
        # get smiles from compound ID
        smi = cache.get('cid_strc')[cmp_id]['smiles']
        # add the compound x sto coeff in the template reaction
        trans_res['left'].extend([smi]*tmpl_rxn['left'][cmp_id])
    # add elements already in input transfo because sto coeff could be > 1 in template reaction
    for cmp_id in detected_compounds['common']['left']:
        # get smiles from compound ID
        smi = cache.get('cid_strc')[cmp_id]['smiles']
        # add the compound x-1 (already in input transfo) sto coeff in the template reaction
        trans_res['left'].extend([smi]*(tmpl_rxn['left'][cmp_id]-1))

    # RIGHT
    for cmp_id in detected_compounds['toadd']['right']:
        # get smiles from compound ID
        smi = cache.get('cid_strc')[cmp_id]['smiles']
        # add the compound x sto coeff in the template reaction
        trans_res['right'].extend([cmp_id]*tmpl_rxn['right'][cmp_id])
    # add elements already in input transfo because sto coeff could be > 1 in template reaction
    for cmp_id in detected_compounds['common']['right']:
        # get smiles from compound ID
        smi = cache.get('cid_strc')[cmp_id]['smiles']
        # add the compound x-1 (already in input transfo) sto coeff in the template reaction
        trans_res['right'].extend([cmp_id]*(tmpl_rxn['right'][cmp_id]-1))

    return trans_res