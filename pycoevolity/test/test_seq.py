#! /usr/bin/env python

import unittest
import os
import sys
import logging

from pycoevolity import seq
from pycoevolity.test import TestLevel
from pycoevolity.test.support.pycoevolity_test_case import PycoevolityTestCase

_LOG = logging.getLogger(__name__)


class GetPropSharedAndDiffTestCase(unittest.TestCase):
    def test_unaligned(self):
        seq_a = "AAGTC-?N-?"
        seq_b = "AGTC-?N-?"
        self.assertRaises(Exception, seq.get_overlap_and_diff, seq_a, seq_b)

    def test_one_missing(self):
        seq_a = "?????"
        seq_b = "AGTCG"
        p_shared, p_diff = seq.get_overlap_and_diff(seq_a, seq_b)
        self.assertTrue(p_shared == 0.0)
        self.assertTrue(p_diff is None)

    def test_both_missing(self):
        seq_a = "?????"
        seq_b = "?????"
        p_shared, p_diff = seq.get_overlap_and_diff(seq_a, seq_b)
        self.assertTrue(p_shared == 0.0)
        self.assertTrue(p_diff is None)

    def test_mutual_missing(self):
        seq_a = "?????ACGT"
        seq_b = "AGTCG????"
        p_shared, p_diff = seq.get_overlap_and_diff(seq_a, seq_b)
        self.assertTrue(p_shared == 0.0)
        self.assertTrue(p_diff is None)

    def test_identical(self):
        seq_a = "ACGT"
        seq_b = "ACGT"
        p_shared, p_diff = seq.get_overlap_and_diff(seq_a, seq_b)
        self.assertTrue(p_shared == 1.0)
        self.assertTrue(p_diff == 0.0)

    def test_fully_shared_diffs(self):
        seq_a = "ACGTG"
        seq_b = "ACGTC"
        p_shared, p_diff = seq.get_overlap_and_diff(seq_a, seq_b)
        self.assertTrue(p_shared == 1.0)
        self.assertTrue(p_diff == (1.0/5.0))

    def test_partial_shared_diffs(self):
        seq_a = "AC?GTG"
        seq_b = "ACGGTC"
        p_shared, p_diff = seq.get_overlap_and_diff(seq_a, seq_b)
        self.assertTrue(p_shared == (5.0/6.0))
        self.assertTrue(p_diff == (1.0/5.0))
        
if __name__ == '__main__':
    unittest.main()

