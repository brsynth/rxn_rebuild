"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from unittest  import TestCase
from rr_cache import rrCache
from rxn_rebuild.rxn_rebuild import (
    rebuild_rxn,
    complete_transfo
)
from brs_utils import create_logger


class Test(TestCase):

    # Set attributes
    logger = create_logger(__name__, 'INFO')
    cache = rrCache(
        attrs=['rr_reactions', 'template_reactions', 'cid_strc']
    )

    def test_all_cmpds_ok(self):
        rule_id = 'RR-02-0250d458c4991a7d-02-F'
        transfo = '[H][O][C](=[O])[C](=[O])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[O][H]>>[H]OC(=O)C(=O)C([H])(O[H])C([H])(O[H])C([H])([H])C([H])=O.[H]O[H].[H]O[H]'
        completed_transfos = rebuild_rxn(
                  cache = self.cache,
            rxn_rule_id = rule_id,
                transfo = transfo,
                 logger = self.logger
        )
        self.assertEqual(
            completed_transfos,
            {
                'MNXR94682': {
                    'full_transfo': {
                        'left': {
                            '[H][O][C](=[O])[C](=[O])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[O][H]': 1.0
                        },
                        'right': {
                            '[H]OC(=O)C(=O)C([H])(O[H])C([H])(O[H])C([H])([H])C([H])=O': 1.0,
                            '[H]O[H]': 2.0
                        }
                    },
                    'added_cmpds': {
                        'left': {},
                        'right': {},
                        'left_nostruct': {},
                        'right_nostruct': {}
                    },
                    'sep_side': '>>',
                    'sep_cmpd': '.'
                }
            }
        )

    def test_cmpds_nostruct(self):
        rule_id = 'RR-02-5594c358cc56f547-02-F'
        transfo = '[H][O][C]([H])([H])[C]([H])([H])[C]([H])([H])[O][H]>>[H+].[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H].[H][O][C]([H])([H])[C]([H])([H])[C]([H])=[O]'
        completed_transfos = rebuild_rxn(
                  cache = self.cache,
            rxn_rule_id = rule_id,
                transfo = transfo,
                 logger = self.logger
        )
        self.assertEqual(
            completed_transfos,
            {
                "MNXR94690": {
                    "full_transfo": {
                        "left": {
                            "[H][O][C]([H])([H])[C]([H])([H])[C]([H])([H])[O][H]": 1.0,
                            "O": 1.0
                        },
                        "right": {
                            "[H+]": 2.0,
                            "[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]": 1.0,
                            "[H][O][C]([H])([H])[C]([H])([H])[C]([H])=[O]": 1.0
                        }
                    },
                    "added_cmpds": {
                        "left": {
                            "MNXM2": {
                                "stoichio": 1.0,
                                "formula": "H2O",
                                "smiles": "O",
                                "inchi": "InChI=1S/H2O/h1H2",
                                "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                                "cid": "MNXM2",
                                "name": "H2O"
                            }
                        },
                        "right": {
                            "MNXM1": {
                                "stoichio": 1.0,
                                "formula": "H",
                                "smiles": "[H+]",
                                "inchi": "InChI=1S/p+1",
                                "inchikey": "GPRLSGONYQIRFK-UHFFFAOYSA-N",
                                "cid": "MNXM1",
                                "name": "H(+)"
                            }
                        },
                        "left_nostruct": {},
                        "right_nostruct": {}
                    },
                    "sep_side": ">>",
                    "sep_cmpd": "."
                },
                "MNXR119958": {
                    "full_transfo": {
                        "left": {
                            "[H][O][C]([H])([H])[C]([H])([H])[C]([H])([H])[O][H]": 1.0,
                            "O": 1.0,
                            "MNXM8975": 1.0
                        },
                        "right": {
                            "[H+]": 1.0,
                            "[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]": 1.0,
                            "[H][O][C]([H])([H])[C]([H])([H])[C]([H])=[O]": 1.0,
                            "MNXM36": 1.0
                        }
                    },
                    "added_cmpds": {
                        "left": {
                            "MNXM2": {
                                "stoichio": 1.0,
                                "formula": "H2O",
                                "smiles": "O",
                                "inchi": "InChI=1S/H2O/h1H2",
                                "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                                "cid": "MNXM2",
                                "name": "H2O"
                            }
                        },
                        "right": {},
                        "left_nostruct": {
                            "MNXM8975": {
                                "stoichio": 1.0,
                                "cid": "MNXM8975"
                            }
                        },
                        "right_nostruct": {
                            "MNXM36": {
                                "stoichio": 1.0,
                                "cid": "MNXM36"
                            }
                        }
                    },
                    "sep_side": ">>",
                    "sep_cmpd": "."
                },
                "MNXR125838": {
                    "full_transfo": {
                        "left": {
                            "[H][O][C]([H])([H])[C]([H])([H])[C]([H])([H])[O][H]": 1.0,
                            "O": 1.0,
                            "MNXM35": 1.0
                        },
                        "right": {
                            "[H+]": 1.0,
                            "[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]": 1.0,
                            "[H][O][C]([H])([H])[C]([H])([H])[C]([H])=[O]": 1.0,
                            "MNXM24": 1.0
                        }
                    },
                    "added_cmpds": {
                        "left": {
                            "MNXM2": {
                                "stoichio": 1.0,
                                "formula": "H2O",
                                "smiles": "O",
                                "inchi": "InChI=1S/H2O/h1H2",
                                "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                                "cid": "MNXM2",
                                "name": "H2O"
                            }
                        },
                        "right": {},
                        "left_nostruct": {
                            "MNXM35": {
                                "stoichio": 1.0,
                                "cid": "MNXM35"
                            }
                        },
                        "right_nostruct": {
                            "MNXM24": {
                                "stoichio": 1.0,
                                "cid": "MNXM24"
                            }
                        }
                    },
                    "sep_side": ">>",
                    "sep_cmpd": "."
                },
                "MNXR128882": {
                    "full_transfo": {
                        "left": {
                            "[H][O][C]([H])([H])[C]([H])([H])[C]([H])([H])[O][H]": 1.0,
                            "O": 1.0,
                            "MNXM35": 1.0
                        },
                        "right": {
                            "[H+]": 1.0,
                            "[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]": 1.0,
                            "[H][O][C]([H])([H])[C]([H])([H])[C]([H])=[O]": 1.0,
                            "MNXM24": 1.0
                        }
                    },
                    "added_cmpds": {
                        "left": {
                            "MNXM2": {
                                "stoichio": 1.0,
                                "formula": "H2O",
                                "smiles": "O",
                                "inchi": "InChI=1S/H2O/h1H2",
                                "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                                "cid": "MNXM2",
                                "name": "H2O"
                            }
                        },
                        "right": {},
                        "left_nostruct": {
                            "MNXM35": {
                                "stoichio": 1.0,
                                "cid": "MNXM35"
                            }
                        },
                        "right_nostruct": {
                            "MNXM24": {
                                "stoichio": 1.0,
                                "cid": "MNXM24"
                            }
                        }
                    },
                    "sep_side": ">>",
                    "sep_cmpd": "."
                }
            }
        )

    def test_all_cmpds_ok_wocache(self):
        rule_id = 'RR-02-0250d458c4991a7d-02-F'
        transfo = '[H][O][C](=[O])[C](=[O])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[O][H]>>[H]OC(=O)C(=O)C([H])(O[H])C([H])(O[H])C([H])([H])C([H])=O.[H]O[H].[H]O[H]'
        completed_transfos = rebuild_rxn(
            rxn_rule_id = rule_id,
                transfo = transfo,
                 logger = self.logger
        )
        self.assertEqual(
            completed_transfos,
            {
                "MNXR94682": {
                    "full_transfo": {
                        "left": {
                            "[H][O][C](=[O])[C](=[O])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[O][H]": 1.0
                        },
                        "right": {
                            "[H]OC(=O)C(=O)C([H])(O[H])C([H])(O[H])C([H])([H])C([H])=O": 1.0,
                            "[H]O[H]": 2.0
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
            }
        )

    def test_forward_direction(self):
        rule_id = 'RR-02-a0cc0be463ff412f-16-F'
        transfo = '[H]Oc1c([H])c([H])c([H])c([H])c1O[H].O=O>>[H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]'
        direction = 'forward'
        self.assertEqual(
            rebuild_rxn(
                rxn_rule_id = rule_id,
                transfo = transfo,
                direction = direction,
                logger = self.logger
            ),
            {
                "MNXR96458": {
                    "full_transfo": {
                        'left': {
                            '[H]Oc1c([H])c([H])c([H])c([H])c1O[H]': 1.0,
                            'O=O': 1.0
                        },
                        'right': {
                            '[H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]': 1.0,
                            '[H+]': 2.0
                        },
                    },
                    'sep_side': '>>',
                    'sep_cmpd': '.',
                    "added_cmpds": {
                        "left": {},
                        "right": {
                            "MNXM1": {
                                "stoichio": 2.0,
                                "formula": "H",
                                "smiles": "[H+]",
                                "inchi": "InChI=1S/p+1",
                                "inchikey": "GPRLSGONYQIRFK-UHFFFAOYSA-N",
                                "cid": "MNXM1",
                                "name": "H(+)"
                            }
                        },
                        "left_nostruct": {},
                        "right_nostruct": {}
                    }
                }
            }
        )

    def test_complete_transfo(self):
        trans_input = {'left': {'[H]OC(=O)C(=O)C([H])([H])[H]': 1, '[H]OO[H]': 1}, 'right': {'[H]OC(=O)C([H])(O[H])C([H])([H])[H]':1}, 'format': 'smiles', 'sep_side': '>>', 'sep_cmpd': '.'}
        rxn_rule = {'rule_id': 'RR-02-0f2f669fd09feecb-12-F', 'rule_score': 0.3819283206707785, 'reac_id': 'MNXR133653', 'subs_id': 'MNXM179', 'rel_direction': 1, 'left': {'MNXM179': 1}, 'right': {'MNXM22': 1, 'MNXM23': 1}}
        tmpl_rxn = {'left': {'MNXM179': 1, 'MNXM4': 1}, 'right': {'MNXM22': 1, 'MNXM23': 1}, 'direction': 1, 'main_left': ['MNXM179'], 'main_right': ['MNXM22', 'MNXM23']}
        tmpl_rxn_id = 'MNXR133653'
        compounds = {'left': {'MNXM4': {'stoichio': 1, 'formula': 'O2', 'smiles': 'O=O', 'inchi': 'InChI=1S/O2/c1-2', 'inchikey': 'MYMOFIZGZYHOMD-UHFFFAOYSA-N', 'cid': 'MNXM4', 'name': 'O2'}}, 'right': {}, 'left_nostruct': {}, 'right_nostruct': {}}
        # FORWARD (x2 to see if results are same)
        for i in range(2):
            compl_transfo = complete_transfo(
                trans_input=trans_input,
                rxn_rule=rxn_rule,
                tmpl_rxn=tmpl_rxn,
                tmpl_rxn_id=tmpl_rxn_id,
                compounds=compounds,
                direction='forward'
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
                            "[H]OC(=O)C([H])(O[H])C([H])([H])[H]": 1,
                            "MNXM4": 1
                        }
                    },
                    "added_cmpds": {
                        "left": {},
                        "right": {},
                        "left_nostruct": {},
                        "right_nostruct": {
                            "MNXM4": {
                                "stoichio": 1,
                                "cid": "MNXM4"
                            }
                        }
                    },
                    "sep_side": ">>",
                    "sep_cmpd": "."
                }
            )
        # REVERSE (x2 to see if results are same)
        for i in range(2):
            compl_transfo = complete_transfo(
                trans_input=trans_input,
                rxn_rule=rxn_rule,
                tmpl_rxn=tmpl_rxn,
                tmpl_rxn_id=tmpl_rxn_id,
                compounds=compounds,
                direction='reverse'
            )
            self.assertDictEqual(
                compl_transfo,
                {
                    "full_transfo": {
                        "left": {
                            "[H]OC(=O)C(=O)C([H])([H])[H]": 1,
                            "[H]OO[H]": 1,
                            "MNXM4": 1
                        },
                        "right": {
                            "[H]OC(=O)C([H])(O[H])C([H])([H])[H]": 1
                        }
                    },
                    "added_cmpds": {
                        "left": {},
                        "right": {},
                        "left_nostruct": {
                            "MNXM4": {
                                "stoichio": 1,
                                "cid": "MNXM4"
                            }
                        },
                        "right_nostruct": {}
                    },
                    "sep_side": ">>",
                    "sep_cmpd": "."
                }
            )

    def test_complete_transfo_with_compds_to_ignore(self):
        trans_input = {'left': {'[H]OC(=O)C(=O)C([H])([H])[H]': 1, '[H]OO[H]': 1}, 'right': {'[H]OC(=O)C([H])(O[H])C([H])([H])[H]':1}, 'format': 'smiles', 'sep_side': '>>', 'sep_cmpd': '.'}
        rxn_rule = {'rule_id': 'RR-02-0f2f669fd09feecb-12-F', 'rule_score': 0.3819283206707785, 'reac_id': 'MNXR133653', 'subs_id': 'MNXM179', 'rel_direction': 1, 'left': {'MNXM179': 1}, 'right': {'MNXM22': 1, 'MNXM23': 1}}
        tmpl_rxn = {'left': {'MNXM179': 1, 'MNXM4': 1}, 'right': {'MNXM22': 1, 'MNXM23': 1}, 'direction': 1, 'main_left': ['MNXM179'], 'main_right': ['MNXM22', 'MNXM23']}
        tmpl_rxn_id = 'MNXR133653'
        compounds = {'left': {'MNXM4': {'stoichio': 1, 'formula': 'O2', 'smiles': 'O=O', 'inchi': 'InChI=1S/O2/c1-2', 'inchikey': 'MYMOFIZGZYHOMD-UHFFFAOYSA-N', 'cid': 'MNXM4', 'name': 'O2'}}, 'right': {}, 'left_nostruct': {}, 'right_nostruct': {}}
        cmpds_to_ignore = ['MNXM4']
        # FORWARD (x2 to see if results are same)
        for i in range(2):
            compl_transfo = complete_transfo(
                trans_input=trans_input,
                rxn_rule=rxn_rule,
                tmpl_rxn=tmpl_rxn,
                tmpl_rxn_id=tmpl_rxn_id,
                compounds=compounds,
                direction='forward',
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
        # REVERSE (x2 to see if results are same)
        for i in range(2):
            compl_transfo = complete_transfo(
                trans_input=trans_input,
                rxn_rule=rxn_rule,
                tmpl_rxn=tmpl_rxn,
                tmpl_rxn_id=tmpl_rxn_id,
                compounds=compounds,
                direction='reverse',
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