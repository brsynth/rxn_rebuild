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
from rdkit import Chem
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
    logger.debug('INPUT TRANSFORMATION:'+str(dumps(trans_input, indent=4)))

    ## REACTION RULE
    orig_rxn_id, rxn_rule = load_rxn_rule(cache, rxn_rule_id)
    logger.debug('REACTION RULE:'+str(dumps(rxn_rule, indent=4)))

    # Check if the number of structures in SMILES trans_input is the same as in the original reaction
    # Just in right part since rules are always mono-component
    if not check_nb_compounds(
        trans_input['right'],
        rxn_rule['right'],
        logger
    ):
        logger.warning(
            '      Number of compounds different in the original reaction and the input transformation'
        )
        logger.warning('         |- INPUT TRANSFORMATION [RIGHT]: ' + str(trans_input['right']))
        logger.warning('         |- ORIGINAL REACTION [RIGHT]: ' + str(rxn_rule['right']))

    ## ORIGINAL REACTION
    orig_rxn = load_orig_rxn(cache, orig_rxn_id, rxn_rule['rel_direction'])
    logger.debug('ORIGINAL REACTION:'+str(dumps(orig_rxn, indent=4)))

    ## ADD MISSING COMPOUNDS INTO FINAL TRANSFORMATION
    trans_res = complete_input_trans(
        cache,
        trans_input,
        orig_rxn,
        rxn_rule)
    logger.debug('COMPLETED TRANSFORMATION:'+str(dumps(trans_res, indent=4)))

    # Check if the number of structures in SMILES trans_res is the same as in the original reaction
    for side in ['left', 'right']:
        if not check_nb_compounds(
            trans_res[side],
            rxn_rule[side],
            logger
        ):
            logger.warning(
                '      Number of compounds different in the original reaction and the input transformation'
            )
            logger.warning('         |- INPUT TRANSFORMATION ['+side+']: ' + str(trans_input[side]))
            logger.warning('         |- ORIGINAL REACTION ['+side+']: ' + str(rxn_rule[side]))


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


def check_nb_compounds(
    cmp_lst_1: List,
    cmp_lst_2: Dict,
    logger: Logger=getLogger(__file__)
) -> bool:
    return len(cmp_lst_1) == sum(cmp_lst_2.values())


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
        ID of the original reaction.
    rxn_rule: Dict
        Reaction Rule looked for.
    """
    cache.load(['rr_reactions'])
    # Take the first reaction from the list of original reactions
    rxn_id = next(iter(cache.get('rr_reactions')[rxn_rule_id]))
    return rxn_id, cache.get('rr_reactions')[rxn_rule_id][rxn_id]


def load_orig_rxn(
    cache: 'rrCache',
    rxn_id: str,
    rr_direction: int,
    logger: Logger=getLogger(__file__)
) -> Dict:
    """
    Seeks for the original reaction of ID rxn_id in the cache.

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
    orig_rxn: Dict
        Original reaction looked for.
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
    orig_rxn: Dict,
    rxn_rule: Dict,
    logger: Logger=getLogger(__file__)
) -> Dict:
    """
    Detect compounds that are different and the same between the original reaction and the reaction rule.

    Parameters
    ----------
    orig_rxn: Dict
        Orignial reaction.
    rxn_rule: Dict
        Reaction rule.
    logger : Logger
        The logger object.

    Returns
    -------
    det_cmp: Dict
        Compounds:
            - different between original reaction and the reaction rule
            - same between original reaction and the reaction rul
    """
    # building the sets
    rxn_left_set = set(orig_rxn['left'].keys())
    rxn_right_set = set(orig_rxn['right'].keys())
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
    orig_rxn: Dict,
    rxn_rule: Dict,
    logger: Logger=getLogger(__file__)
) -> Dict:
    """
    Complete input transformation from original reaction and reaction rule.

    Parameters
    ----------
    cache: rrCache
        Pre-computed data.
    trans_input: Dict
        Transformation to complete
    orig_rxn: Dict
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

    ## DETECT IDENTICAL COMPOUNDS IN BOTH REACTION RULE AND ORIGINAL REACTION
    detected_compounds = detect_compounds(orig_rxn, rxn_rule)

    # LEFT
    for cmp_id in detected_compounds['toadd']['left']:
        # get smiles from compound ID
        smi = cache.get('cid_strc')[cmp_id]['smiles']
        # add the compound x sto coeff in the original reaction
        trans_res['left'].extend([smi]*orig_rxn['left'][cmp_id])
    # add elements already in input transfo because sto coeff could be > 1 in original reaction
    for cmp_id in detected_compounds['common']['left']:
        # get smiles from compound ID
        smi = cache.get('cid_strc')[cmp_id]['smiles']
        # add the compound x-1 (already in input transfo) sto coeff in the original reaction
        trans_res['left'].extend([smi]*(orig_rxn['left'][cmp_id]-1))

    # RIGHT
    for cmp_id in detected_compounds['toadd']['right']:
        # get smiles from compound ID
        smi = cache.get('cid_strc')[cmp_id]['smiles']
        # add the compound x sto coeff in the original reaction
        trans_res['right'].extend([cmp_id]*orig_rxn['right'][cmp_id])
    # add elements already in input transfo because sto coeff could be > 1 in original reaction
    for cmp_id in detected_compounds['common']['right']:
        # get smiles from compound ID
        smi = cache.get('cid_strc')[cmp_id]['smiles']
        # add the compound x-1 (already in input transfo) sto coeff in the original reaction
        trans_res['right'].extend([cmp_id]*(orig_rxn['right'][cmp_id]-1))

    return trans_res