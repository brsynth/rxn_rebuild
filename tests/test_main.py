"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from unittest  import TestCase
from os import path as os_path
from json import load as json_load
from rr_cache import rrCache
from rxn_rebuild.rxn_rebuild import (
    rebuild_rxn,
    complete_transfo
)
from brs_utils import create_logger


HERE = os_path.dirname(os_path.abspath(__file__))
DATA_PATH = os_path.join(HERE, 'data')

class Test(TestCase):

    # Set attributes
    logger = create_logger(__name__, 'INFO')
    cache = rrCache(
        # attrs=['rr_reactions', 'template_reactions', 'cid_strc']
        cspace='mnx3.1',
        interactive=False,
        logger=logger
    )
    # cache_rr2026 = 

    def test_all_cmpds_ok(self):
        rule_id = 'RR-02-0250d458c4991a7d-02-F'
        transfo = '[H][O][C](=[O])[C](=[O])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[O][H]>>[H]OC(=O)C(=O)C([H])(O[H])C([H])(O[H])C([H])([H])C([H])=O.[H]O[H].[H]O[H]'
        completed_transfos = rebuild_rxn(
                cache = self.cache,
                rxn_rule_id = rule_id,
                transfo = transfo,
                cspace = 'mnx4.4',
                cspace_type = 'legacy',
                logger = self.logger
        )
        self.assertEqual(
            completed_transfos,
            json_load(open(os_path.join(DATA_PATH, 'all_cmpds_ok.json'), 'r'))
        )

    def test_cmpds_nostruct(self):
        rule_id = 'RR-02-5594c358cc56f547-02-F'
        transfo = '[H][O][C]([H])([H])[C]([H])([H])[C]([H])([H])[O][H]>>[H+].[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H].[H][O][C]([H])([H])[C]([H])([H])[C]([H])=[O]'
        completed_transfos = rebuild_rxn(
                  cache = self.cache,
            rxn_rule_id = rule_id,
                transfo = transfo,
                cspace = 'mnx4.4',
                cspace_type = 'legacy',
                logger = self.logger
        )
        self.assertEqual(
            completed_transfos,
            json_load(open(os_path.join(DATA_PATH, 'cmpds_nostruct.json'), 'r'))
        )

    def test_all_cmpds_ok_wocache(self):
        rule_id = 'RR-02-0250d458c4991a7d-02-F'
        transfo = '[H][O][C](=[O])[C](=[O])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[O][H]>>[H]OC(=O)C(=O)C([H])(O[H])C([H])(O[H])C([H])([H])C([H])=O.[H]O[H].[H]O[H]'
        completed_transfos = rebuild_rxn(
            rxn_rule_id = rule_id,
                transfo = transfo,
                cspace = 'mnx4.4',
                cspace_type = 'legacy',
                logger = self.logger
        )
        self.assertEqual(
            completed_transfos,
            json_load(open(os_path.join(DATA_PATH, 'all_cmpds_ok_wocache.json'), 'r'))
        )

    def test_complete_transfo(self):
        trans_input = {
            'left': {
                'Cc1ccc(C(C)C)cc1O': 1
            },
            'right': {
                'CC1=CCC(C(C)C)=CC1': 1,
                'O=O': 1
            },
            'format': 'smiles',
            'sep_side': '>>',
            'sep_cmpd': '.'
        }
        rxn_rule = {
            "left": {
                "CHEBI:3440": 1
            },
            "left_excluded": [
                "CHEBI:15377",
                "CHEBI:15377",
                "CHEBI:15377",
                "CHEBI:15378",
                "CHEBI:15378",
                "CHEBI:58210",
                "CHEBI:58210"
            ],
            "reac_id": "RHEA:67404",
            "rel_direction": -1,
            "right": {
                "CHEBI:10577": 1,
                "CHEBI:15379": 1
            },
            "right_excluded": [
                "CHEBI:15379",
                "CHEBI:57618",
                "CHEBI:57618"
            ],
            "rule_id": "RR:03-6A67ED-190DBF-97D570",
            "rule_score": 1.0,
            "subs_id": "CHEBI:3440"
        }
        tmpl_rxn = {
            "direction": 0,
            "left": {
                "CHEBI:10577": 1,
                "CHEBI:15379": 2,
                "CHEBI:57618": 2
            },
            "main_left": "CHEBI:10577",
            "main_right": "CHEBI:3440",
            "right": {
                "CHEBI:15377": 3,
                "CHEBI:15378": 2,
                "CHEBI:3440": 1,
                "CHEBI:58210": 2
            }
        }
        tmpl_rxn_id = 'RHEA:67404'
        compounds = {}
        compl_transfo = complete_transfo(
            trans_input=trans_input,
            rxn_rule=rxn_rule,
            tmpl_rxn=tmpl_rxn,
            tmpl_rxn_id=tmpl_rxn_id,
            compounds=compounds,
            legacy=False
        )
        self.assertDictEqual(
            compl_transfo,
            {
                'full_transfo': {
                    'right': {
                        'CC1=CCC(C(C)C)=CC1': 1.0,
                        'O=O': 1.0,
                        'CHEBI:15379': 1,
                        'CHEBI:57618': 2
                    },
                    'left': {
                        'Cc1ccc(C(C)C)cc1O': 1.0,
                        'CHEBI:15377': 3,
                        'CHEBI:15378': 2,
                        'CHEBI:58210': 2
                    }
                },
                'added_cmpds': {
                    'right': {
                        'CHEBI:15379': 1,
                        'CHEBI:57618': 2
                    },
                    'left': {
                        'CHEBI:15377': 3,
                        'CHEBI:15378': 2,
                        'CHEBI:58210': 2
                    }
                },
                'sep_side': '>>',
                'sep_cmpd': '.'
            }
        )

    def test_complete_transfo_with_compds_to_ignore(self):
        trans_input = {'left': {'[H]OC(=O)C(=O)C([H])([H])[H]': 1, '[H]OO[H]': 1}, 'right': {'[H]OC(=O)C([H])(O[H])C([H])([H])[H]':1}, 'format': 'smiles', 'sep_side': '>>', 'sep_cmpd': '.'}
        rxn_rule = {'rule_id': 'RR-02-0f2f669fd09feecb-12-F', 'rule_score': 0.3819283206707785, 'reac_id': 'MNXR133653', 'subs_id': 'MNXM179', 'rel_direction': 1, 'left': {'MNXM179': 1}, 'right': {'MNXM22': 1, 'MNXM23': 1}}
        tmpl_rxn = {'left': {'MNXM179': 1, 'MNXM4': 1}, 'right': {'MNXM22': 1, 'MNXM23': 1}, 'direction': 1, 'main_left': ['MNXM179'], 'main_right': ['MNXM22', 'MNXM23']}
        tmpl_rxn_id = 'MNXR133653'
        compounds = {'left': {'MNXM4': {'stoichio': 1, 'formula': 'O2', 'smiles': 'O=O', 'inchi': 'InChI=1S/O2/c1-2', 'inchikey': 'MYMOFIZGZYHOMD-UHFFFAOYSA-N', 'cid': 'MNXM4', 'name': 'O2'}}, 'right': {}, 'left_nostruct': {}, 'right_nostruct': {}}
        cmpds_to_ignore = ['MNXM4']
        compl_transfo = complete_transfo(
            trans_input=trans_input,
            rxn_rule=rxn_rule,
            tmpl_rxn=tmpl_rxn,
            tmpl_rxn_id=tmpl_rxn_id,
            compounds=compounds,
            legacy=True,
            cmpds_to_ignore=cmpds_to_ignore
        )
        self.assertDictEqual(
            compl_transfo,
            {
                "full_transfo": {
                    "left": {
                        "[H]OC(=O)C(=O)C([H])([H])[H]": 1,
                        "[H]OO[H]": 1
                    },
                    "right": {
                        "[H]OC(=O)C([H])(O[H])C([H])([H])[H]": 1
                    }
                },
                "added_cmpds": {
                    "left": {},
                    "right": {},
                    "left_nostruct": {},
                    "right_nostruct": {}
                },
                "sep_side": ">>",
                "sep_cmpd": "."
            }
        )
