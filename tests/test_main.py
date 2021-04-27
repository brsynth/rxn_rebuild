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

    def test_1(self):
        rule_id = 'RR-02-0250d458c4991a7d-02-F'
        transfo = '[H][O][C](=[O])[C](=[O])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[O][H]>>[H]OC(=O)C(=O)C([H])(O[H])C([H])(O[H])C([H])([H])C([H])=O.[H]O[H].[H]O[H]'
        compl_rxn_smi = rebuild_rxn(
              cache = self.cache,
        rxn_rule_id = rule_id,
          trans_smi = transfo,
             logger = self.logger
        )
        self.assertEqual(
            compl_rxn_smi,
            '[H][O][C](=[O])[C](=[O])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[O][H]>>[H]OC(=O)C(=O)C([H])(O[H])C([H])(O[H])C([H])([H])C([H])=O.[H]O[H].[H]O[H]'
            )
