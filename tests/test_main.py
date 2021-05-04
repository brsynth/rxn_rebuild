"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from unittest  import TestCase
from rr_cache import rrCache
from rxn_rebuild.rxn_rebuild import (
    rebuild_rxn,
)
from brs_utils import create_logger


class Test(TestCase):

    # Set attributes
    logger = create_logger(__name__, 'DEBUG')
    cache = rrCache(
        db='file',
        attrs=None
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
                "MNXR94682": {
                    "added_cmpds": {
                        "left": {},
                        "right": {},
                        "left_nostruct": {},
                        "right_nostruct": {}
                    },
                    "full_transfo": "[H][O][C](=[O])[C](=[O])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[O][H]>>[H]OC(=O)C(=O)C([H])(O[H])C([H])(O[H])C([H])([H])C([H])=O.[H]O[H].[H]O[H]"
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
                    "added_cmpds": {
                        "left": {
                            "MNXM2": {
                                "stoichio": 1,
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
                                "stoichio": 1,
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
                    "full_transfo": "[H][O][C]([H])([H])[C]([H])([H])[C]([H])([H])[O][H].O>>[H+].[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H].[H][O][C]([H])([H])[C]([H])([H])[C]([H])=[O].[H+]"
                },
                "MNXR119958": {
                    "added_cmpds": {
                        "left": {
                            "MNXM2": {
                                "stoichio": 1,
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
                                "stoichio": 1,
                                "cid": "MNXM8975"
                            }
                        },
                        "right_nostruct": {
                            "MNXM36": {
                                "stoichio": 1,
                                "cid": "MNXM36"
                            }
                        }
                    }
                },
                "MNXR125838": {
                    "added_cmpds": {
                        "left": {
                            "MNXM2": {
                                "stoichio": 1,
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
                                "stoichio": 1,
                                "cid": "MNXM35"
                            }
                        },
                        "right_nostruct": {
                            "MNXM24": {
                                "stoichio": 1,
                                "cid": "MNXM24"
                            }
                        }
                    }
                },
                "MNXR128882": {
                    "added_cmpds": {
                        "left": {
                            "MNXM2": {
                                "stoichio": 1,
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
                                "stoichio": 1,
                                "cid": "MNXM35"
                            }
                        },
                        "right_nostruct": {
                            "MNXM24": {
                                "stoichio": 1,
                                "cid": "MNXM24"
                            }
                        }
                    }
                }
            }
        )
