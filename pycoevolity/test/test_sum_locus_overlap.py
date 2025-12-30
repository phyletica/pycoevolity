#! /usr/bin/env python

import unittest
import os
import sys
import logging

from pycoevolity.cli import sum_locus_overlap 
from pycoevolity.test import TestLevel
from pycoevolity.test import test_utils 
from pycoevolity.test.support import package_paths
from pycoevolity.test.support.pycoevolity_test_case import PycoevolityTestCase

_LOG = logging.getLogger(__name__)


class SumLocusOverlapTestCase(PycoevolityTestCase):
    def setUp(self):
        self.set_up()

    def tearDown(self):
        self.tear_down()

    def test_with_json(self):
        ipyrad_loci_path = package_paths.data_path('all-sals.loci')
        ipyrad_json_path = package_paths.data_path('all-sals.json')
        test_dir = self.get_test_subdir(prefix='temp-sum-locus-overlap')
        test_prefix = os.path.join(test_dir, self.test_id) + "-"

        min_overlap = 0.2

        args = [
            "--json-path", ipyrad_json_path,
            "--min-overlap", f"{min_overlap}",
            "--prefix", test_prefix,
            ipyrad_loci_path,
        ]

        sum_locus_overlap.main(argv = args)

        pairwise_results_file = test_prefix + "pycoevolity-pairwise-overlap-and-div.tsv"
        sample_results_file = test_prefix + "pycoevolity-sample-mean-overlap.tsv"

        pairwise_lines = []
        with open(pairwise_results_file, "r") as in_stream:
            for line in in_stream:
                pairwise_lines.append(line)

        sample_lines = []
        with open(sample_results_file, "r") as in_stream:
            for line in in_stream:
                sample_lines.append(line)

        expected_pwise_lines = test_utils.get_expected_pairwise_overlap_div_table(
            ipyrad_loci_path = ipyrad_loci_path,
            min_overlap = min_overlap,
            ipyrad_json_path = ipyrad_json_path,
            stat_key = "reads_passed_filter",
        )

        self.assertEqual(sorted(pairwise_lines), sorted(expected_pwise_lines))

        expected_sample_lines = test_utils.get_expected_per_sample_overlap_table(
            ipyrad_loci_path = ipyrad_loci_path,
            min_overlap = min_overlap,
            ipyrad_json_path = ipyrad_json_path,
            stat_key = "reads_passed_filter",
        )

        self.assertEqual(sorted(sample_lines), sorted(expected_sample_lines))

        self.register_file_system()

if __name__ == '__main__':
    unittest.main()

