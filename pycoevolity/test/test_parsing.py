#! /usr/bin/env python

import unittest
import os
import sys
import math
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import logging

from pycoevolity import parsing
from pycoevolity.test.support.pycoevolity_test_case import PycoevolityTestCase
from pycoevolity.test.support import package_paths

_LOG = logging.getLogger(__name__)

class GetDictFromSpreadsheetTestCase(PycoevolityTestCase):

    def setUp(self):
        self.set_up()

    def tearDown(self):
        self.tear_down()

    def test_tab_with_head(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1_path = self.get_test_path()
        s2_path = self.get_test_path()
        with open(s1_path, 'w') as s1:
            s1.write('{0}\n'.format(sep.join(header)))
            for i in range(3):
                s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        with open(s2_path, 'w') as s2:
            s2.write('{0}\n'.format(sep.join(header)))
            for i in range(3, len(list(d.values())[0])):
                s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        ret = parsing.get_dict_from_spreadsheets([s1_path, s2_path], sep = sep)
        self.assertEqual(d, ret)

    def test_tab_without_head(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1_path = self.get_test_path()
        s2_path = self.get_test_path()
        with open(s1_path, 'w') as s1:
            for i in range(3):
                s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        with open(s2_path, 'w') as s2:
            for i in range(3, len(list(d.values())[0])):
                s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        ret = parsing.get_dict_from_spreadsheets([s1_path, s2_path], sep = sep, header = header)
        self.assertEqual(d, ret)

    def test_comma_with_head(self):
        sep = ','
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1_path = self.get_test_path()
        s2_path = self.get_test_path()
        with open(s1_path, 'w') as s1:
            s1.write('{0}\n'.format(sep.join(header)))
            for i in range(3):
                s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        with open(s2_path, 'w') as s2:
            s2.write('{0}\n'.format(sep.join(header)))
            for i in range(3, len(list(d.values())[0])):
                s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        ret = parsing.get_dict_from_spreadsheets([s1_path, s2_path], sep = sep)
        self.assertEqual(d, ret)

    def test_comma_without_head(self):
        sep = ','
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1_path = self.get_test_path()
        s2_path = self.get_test_path()
        with open(s1_path, 'w') as s1:
            for i in range(3):
                s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        with open(s2_path, 'w') as s2:
            for i in range(3, len(list(d.values())[0])):
                s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        ret = parsing.get_dict_from_spreadsheets([s1_path, s2_path], sep = sep, header = header)
        self.assertEqual(d, ret)

    def test_error_mismatch_headers(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        header2 = ['PRI_Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1_path = self.get_test_path()
        s2_path = self.get_test_path()
        with open(s1_path, 'w') as s1:
            s1.write('{0}\n'.format(sep.join(header)))
            for i in range(3):
                s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        with open(s2_path, 'w') as s2:
            s2.write('{0}\n'.format(sep.join(header2)))
            for i in range(3, len(list(d.values())[0])):
                s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        self.assertRaises(Exception, parsing.get_dict_from_spreadsheets,
                [s1_path, s2_path], sep = sep)

    def test_error_missing_header(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s1_path = self.get_test_path()
        s2_path = self.get_test_path()
        with open(s1_path, 'w') as s1:
            for i in range(3):
                s1.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        with open(s2_path, 'w') as s2:
            for i in range(3, len(list(d.values())[0])):
                s2.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        self.assertRaises(Exception, parsing.get_dict_from_spreadsheets,
                [s1_path, s2_path], sep = sep)

class DictLineIterTestCase(unittest.TestCase):

    def test_tab_with_head(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s = StringIO()
        s.write('{0}\n'.format(sep.join(header)))
        for i in range(len(list(d.values())[0])):
            s.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        s.seek(0)
        r = StringIO()
        for row in parsing.dict_line_iter(d, sep = sep, header = header):
            r.write(row)
        self.assertEqual(s.getvalue(), r.getvalue())

    def test_tab_without_head(self):
        sep = '\t'
        header = ['PRI.Psi', 'pi.net.1', 'pi.net.2']
        header = sorted(header)
        d = dict(zip(header, [['1','5','3', '2', '6'],
                ['0.2', '0.12', '0.11', '0.33', '0.29'],
                ['0.001', '0.0043', '0.0002', '0.0', '0.0036']]))
        s = StringIO()
        s.write('{0}\n'.format(sep.join(header)))
        for i in range(len(list(d.values())[0])):
            s.write('{0}\n'.format(sep.join([d[h][i] for h in header])))
        s.seek(0)
        r = StringIO()
        for row in parsing.dict_line_iter(d, sep = sep):
            r.write(row)
        self.assertEqual(s.getvalue(), r.getvalue())

class EcoevolityStdOutTestCase(PycoevolityTestCase):
    def setUp(self):
        self.set_up()
        self.restart_std_out_path = package_paths.data_path(
                'ecoevolity-std-out.txt')

    def tearDown(self):
        self.tear_down()

    def test_parsing_restart_std_out(self):
        t = parsing.EcoevolityStdOut(self.restart_std_out_path)
        self.assertEqual(t.number_of_comparisons, 3)
        self.assertEqual(t.run_time, 4812)

        self.assertEqual(t.numbers_of_sites, (500000, 400000, 300000))
        self.assertEqual(t.get_number_of_sites(0), 500000)
        self.assertEqual(t.get_number_of_sites(1), 400000)
        self.assertEqual(t.get_number_of_sites(2), 300000)
        self.assertEqual(t.get_min_number_of_sites(), 300000)
        self.assertEqual(t.get_max_number_of_sites(), 500000)
        self.assertEqual(t.get_mean_number_of_sites(), 400000)

        self.assertEqual(t.numbers_of_variable_sites, (13154, 20005, 29310))
        self.assertEqual(t.get_number_of_variable_sites(0), 13154)
        self.assertEqual(t.get_number_of_variable_sites(1), 20005)
        self.assertEqual(t.get_number_of_variable_sites(2), 29310)
        self.assertEqual(t.get_min_number_of_variable_sites(), 13154)
        self.assertEqual(t.get_max_number_of_variable_sites(), 29310)
        self.assertApproxEqual(t.get_mean_number_of_variable_sites(), 20823.0)

        self.assertEqual(t.numbers_of_patterns, (61, 61, 61))
        self.assertEqual(t.get_number_of_patterns(0), 61)
        self.assertEqual(t.get_number_of_patterns(1), 61)
        self.assertEqual(t.get_number_of_patterns(2), 61)
        self.assertEqual(t.get_min_number_of_patterns(), 61)
        self.assertEqual(t.get_max_number_of_patterns(), 61)
        self.assertApproxEqual(t.get_mean_number_of_patterns(), 61.0)


class MsbayesConfigTestCase(PycoevolityTestCase):

    def setUp(self):
        self.set_up()
        self.cfg_path = package_paths.data_path(
                'msbayes.cfg')
        self.dpp_cfg_path = package_paths.data_path(
                'dpp-msbayes.cfg')
        self.dpp_sps_cfg_path = package_paths.data_path(
                'dpp-msbayes-subs-per-site.cfg')
        self.expected_alignments = {
            'species1': {
                'mt': (
                    (
                        (
                            'species1_mt_org1_pop1',
                            tuple('CTAGGAATTT')
                        ),
                        (
                            'species1_mt_org2_pop1',
                            tuple('CTAGGAATTT')
                        ),
                    ),
                    (
                        (
                            'species1_mt_org3_pop2',
                            tuple('CTAGGAATTC')
                        ),
                        (
                            'species1_mt_org4_pop2',
                            tuple('CTAGGAATTC')
                        ),
                        (
                            'species1_mt_org5_pop2',
                            tuple('CTAGGAATYC')
                        ),
                    ),
                ),
                'nuc': (
                    (
                        (
                            'species1_nuc_org1a_pop1',
                            tuple('AAGTTCTATG')
                        ),
                        (
                            'species1_nuc_org1b_pop1',
                            tuple('AAGTTCTATG')
                        ),
                        (
                            'species1_nuc_org2a_pop1',
                            tuple('AAGTTCTATG')
                        ),
                        (
                            'species1_nuc_org2b_pop1',
                            tuple('AAGTTCTATG')
                        ),
                    ),
                    (
                        (
                            'species1_nuc_org3a_pop2',
                            tuple('GAGTTCTATG')
                        ),
                        (
                            'species1_nuc_org3b_pop2',
                            tuple('GAGTTCTATG')
                        ),
                        (
                            'species1_nuc_org4a_pop2',
                            tuple('GAGTTCTATG')
                        ),
                        (
                            'species1_nuc_org4b_pop2',
                            tuple('GAGTTCTATG')
                        ),
                        (
                            'species1_nuc_org5a_pop2',
                            tuple('GAGTTCTATG')
                        ),
                        (
                            'species1_nuc_org5b_pop2',
                            tuple('GAGTTCTATG')
                        ),
                    ),
                ),
            },
            'species2': {
                'mt': (
                    (
                        (
                            'species2_mt_org1_pop1',
                            tuple('TTAAGAATTT')
                        ),
                        (
                            'species2_mt_org2_pop1',
                            tuple('TTAAGAATTT')
                        ),
                        (
                            'species2_mt_org3_pop1',
                            tuple('TTAAGAATTC')
                        ),
                    ),
                    (
                        (
                            'species2_mt_org4_pop2',
                            tuple('TTAAGAATTC')
                        ),
                        (
                            'species2_mt_org5_pop2',
                            tuple('TTAARAATCC')
                        ),
                    ),
                ),
                'nuc': (
                    (
                        (
                            'species2_nuc_org1a_pop1',
                            tuple('AAGTTCCATA')
                        ),
                        (
                            'species2_nuc_org1b_pop1',
                            tuple('AAGTTCCATA')
                        ),
                        (
                            'species2_nuc_org2a_pop1',
                            tuple('AAGTTCCATA')
                        ),
                        (
                            'species2_nuc_org2b_pop1',
                            tuple('AAGTTCCATA')
                        ),
                        (
                            'species2_nuc_org3a_pop1',
                            tuple('GAGTTCCATA')
                        ),
                        (
                            'species2_nuc_org3b_pop1',
                            tuple('GAGTTCCATA')
                        ),
                    ),
                    (
                        (
                            'species2_nuc_org4a_pop2',
                            tuple('GAGTTCCATA')
                        ),
                        (
                            'species2_nuc_org4b_pop2',
                            tuple('GAGTTCCATA')
                        ),
                        (
                            'species2_nuc_org5a_pop2',
                            tuple('GAGTTCCATA')
                        ),
                        (
                            'species2_nuc_org5b_pop2',
                            tuple('GAGTTCCATA')
                        ),
                    ),
                ),
            },
        }
        self.expected_recoded_alignments = {}
        for taxon, locus_data in self.expected_alignments.items():
            if not taxon in self.expected_recoded_alignments:
                self.expected_recoded_alignments[taxon] = {}
            for locus, align_data in locus_data.items():
                new_seqs = []
                for pop_seqs in align_data:
                    new_pop_seqs = []
                    for seq_label, seq in pop_seqs:
                        new_seq = ''.join(seq)
                        new_seq = new_seq.replace('Y', '?')
                        new_seq = new_seq.replace('R', '?')
                        new_pop_seqs.append(
                            (
                                seq_label,
                                tuple(new_seq),
                            )
                        )
                    new_seqs.append(tuple(new_pop_seqs))
                self.expected_recoded_alignments[taxon][locus] = tuple(new_seqs)
        self.expected_concat_alignments = {
            'species1': [
                        [
                            'species1_mt_org1_pop1_comparison1pop1',
                            'CTAGGAATTT??????????'
                        ],
                        [
                            'species1_mt_org2_pop1_comparison1pop1',
                            'CTAGGAATTT??????????'
                        ],
                        [
                            'species1_mt_org3_pop2_comparison1pop2',
                            'CTAGGAATTC??????????'
                        ],
                        [
                            'species1_mt_org4_pop2_comparison1pop2',
                            'CTAGGAATTC??????????'
                        ],
                        [
                            'species1_mt_org5_pop2_comparison1pop2',
                            'CTAGGAATYC??????????'
                        ],
                        [
                            'species1_nuc_org1a_pop1_comparison1pop1',
                            '??????????AAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org1b_pop1_comparison1pop1',
                            '??????????AAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org2a_pop1_comparison1pop1',
                            '??????????AAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org2b_pop1_comparison1pop1',
                            '??????????AAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org3a_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org3b_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org4a_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org4b_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org5a_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org5b_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
            ],
            'species2': [
                        [
                            'species2_mt_org1_pop1_comparison2pop1',
                            'TTAAGAATTT??????????'
                        ],
                        [
                            'species2_mt_org2_pop1_comparison2pop1',
                            'TTAAGAATTT??????????'
                        ],
                        [
                            'species2_mt_org3_pop1_comparison2pop1',
                            'TTAAGAATTC??????????'
                        ],
                        [
                            'species2_mt_org4_pop2_comparison2pop2',
                            'TTAAGAATTC??????????'
                        ],
                        [
                            'species2_mt_org5_pop2_comparison2pop2',
                            'TTAARAATCC??????????'
                        ],
                        [
                            'species2_nuc_org1a_pop1_comparison2pop1',
                            '??????????AAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org1b_pop1_comparison2pop1',
                            '??????????AAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org2a_pop1_comparison2pop1',
                            '??????????AAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org2b_pop1_comparison2pop1',
                            '??????????AAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org3a_pop1_comparison2pop1',
                            '??????????GAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org3b_pop1_comparison2pop1',
                            '??????????GAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org4a_pop2_comparison2pop2',
                            '??????????GAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org4b_pop2_comparison2pop2',
                            '??????????GAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org5a_pop2_comparison2pop2',
                            '??????????GAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org5b_pop2_comparison2pop2',
                            '??????????GAGTTCCATA'
                        ],
            ],
        }
        self.expected_concat_alignments_recode = {
            'species1': [
                        [
                            'species1_mt_org1_pop1_comparison1pop1',
                            'CTAGGAATTT??????????'
                        ],
                        [
                            'species1_mt_org2_pop1_comparison1pop1',
                            'CTAGGAATTT??????????'
                        ],
                        [
                            'species1_mt_org3_pop2_comparison1pop2',
                            'CTAGGAATTC??????????'
                        ],
                        [
                            'species1_mt_org4_pop2_comparison1pop2',
                            'CTAGGAATTC??????????'
                        ],
                        [
                            'species1_mt_org5_pop2_comparison1pop2',
                            'CTAGGAAT?C??????????'
                        ],
                        [
                            'species1_nuc_org1a_pop1_comparison1pop1',
                            '??????????AAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org1b_pop1_comparison1pop1',
                            '??????????AAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org2a_pop1_comparison1pop1',
                            '??????????AAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org2b_pop1_comparison1pop1',
                            '??????????AAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org3a_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org3b_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org4a_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org4b_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org5a_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
                        [
                            'species1_nuc_org5b_pop2_comparison1pop2',
                            '??????????GAGTTCTATG'
                        ],
            ],
            'species2': [
                        [
                            'species2_mt_org1_pop1_comparison2pop1',
                            'TTAAGAATTT??????????'
                        ],
                        [
                            'species2_mt_org2_pop1_comparison2pop1',
                            'TTAAGAATTT??????????'
                        ],
                        [
                            'species2_mt_org3_pop1_comparison2pop1',
                            'TTAAGAATTC??????????'
                        ],
                        [
                            'species2_mt_org4_pop2_comparison2pop2',
                            'TTAAGAATTC??????????'
                        ],
                        [
                            'species2_mt_org5_pop2_comparison2pop2',
                            'TTAA?AATCC??????????'
                        ],
                        [
                            'species2_nuc_org1a_pop1_comparison2pop1',
                            '??????????AAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org1b_pop1_comparison2pop1',
                            '??????????AAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org2a_pop1_comparison2pop1',
                            '??????????AAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org2b_pop1_comparison2pop1',
                            '??????????AAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org3a_pop1_comparison2pop1',
                            '??????????GAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org3b_pop1_comparison2pop1',
                            '??????????GAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org4a_pop2_comparison2pop2',
                            '??????????GAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org4b_pop2_comparison2pop2',
                            '??????????GAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org5a_pop2_comparison2pop2',
                            '??????????GAGTTCCATA'
                        ],
                        [
                            'species2_nuc_org5b_pop2_comparison2pop2',
                            '??????????GAGTTCCATA'
                        ],
            ],
        }

    def tearDown(self):
        self.tear_down()

    def test_invalid_config(self):
        invalid_cfg_path = package_paths.data_path(
            'ecoevolity-std-out.txt')

        self.assertRaises(Exception, parsing.MsbayesConfig, invalid_cfg_path)

    def test_parsing_msbayes_config(self):
        cfg = parsing.MsbayesConfig(
            self.cfg_path,
            recode_ambig_states_as_missing = False,
        )
        self.assertEqual(cfg.npairs, 2)

        expected_taxa = ["species1", "species2"]
        taxa = sorted(cfg.taxa)
        self.assertEqual(expected_taxa, taxa)

        self.assertEqual(cfg.using_dpp_model(), False)
        self.assertEqual(cfg.using_time_in_subs_per_site(), False)

        expected_eco_model = (
            'global_comparison_settings:\n'
            '    ploidy: 1\n'
            '    genotypes_are_diploid: false\n'
            '    markers_are_dominant: false\n'
            '    constant_sites_removed: false\n'
            '    population_name_delimiter: \" \"\n'
            '    population_name_is_prefix: false\n'
            '    equal_population_sizes: false\n'
            '    parameters:\n'
            '        freq_1:\n'
            '            value: 0.5\n'
            '            estimate: false\n'
            '        mutation_rate:\n'
            '            value: 1.0\n'
            '            estimate: false\n'
            '        root_relative_population_size:\n'
            '            value: 1.0\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 100.0\n'
            '                    scale: 0.01\n'
            '        population_size:\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 5.0\n'
            '                    scale: 0.0002\n'
        )

        s = StringIO()
        cfg.write_ecoevolity_model_settings(s)
        self.assertEqual(expected_eco_model, s.getvalue())

        a = cfg.sample_table.alignments
        alignments = {}
        for taxon, locus_data in a.items():
            if not taxon in alignments:
                alignments[taxon] = {}
            for locus, msb_alignment in locus_data.items():
                alignments[taxon][locus] =  msb_alignment.sequences
        self.assertEqual(alignments, self.expected_alignments)

        concat_alignments = {}
        for comparison_label, data in parsing.Loci.iter_from_msbayes_config(
            cfg,
            remove_triallelic_sites = False,
            convert_to_binary = False,
        ):
            s = StringIO()
            data.write_phylip(s)
            s.seek(0)
            sequences = []
            locus_2_idx = 0
            for i, line in enumerate(s):
                if i < 1:
                    continue
                line = line.strip()
                if line == '':
                    continue
                data = line.split()
                if len(data) == 2:
                    sequences.append(data)
                elif len(data) == 1:
                    sequences[locus_2_idx][1] += data[0]
                    locus_2_idx += 1
                else:
                    raise Exception('Unexpected phylip file line in test')
            concat_alignments[comparison_label] = sequences
        self.assertEqual(self.expected_concat_alignments, concat_alignments)

    def test_parsing_msbayes_config_recode(self):
        cfg = parsing.MsbayesConfig(
            self.cfg_path,
            recode_ambig_states_as_missing = True,
        )
        self.assertEqual(cfg.npairs, 2)

        expected_taxa = ["species1", "species2"]
        taxa = sorted(cfg.taxa)
        self.assertEqual(expected_taxa, taxa)

        self.assertEqual(cfg.using_dpp_model(), False)
        self.assertEqual(cfg.using_time_in_subs_per_site(), False)

        expected_eco_model = (
            'global_comparison_settings:\n'
            '    ploidy: 1\n'
            '    genotypes_are_diploid: false\n'
            '    markers_are_dominant: false\n'
            '    constant_sites_removed: false\n'
            '    population_name_delimiter: \" \"\n'
            '    population_name_is_prefix: false\n'
            '    equal_population_sizes: false\n'
            '    parameters:\n'
            '        freq_1:\n'
            '            value: 0.5\n'
            '            estimate: false\n'
            '        mutation_rate:\n'
            '            value: 1.0\n'
            '            estimate: false\n'
            '        root_relative_population_size:\n'
            '            value: 1.0\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 100.0\n'
            '                    scale: 0.01\n'
            '        population_size:\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 5.0\n'
            '                    scale: 0.0002\n'
        )

        s = StringIO()
        cfg.write_ecoevolity_model_settings(s)
        self.assertEqual(expected_eco_model, s.getvalue())

        a = cfg.sample_table.alignments
        alignments = {}
        for taxon, locus_data in a.items():
            if not taxon in alignments:
                alignments[taxon] = {}
            for locus, msb_alignment in locus_data.items():
                alignments[taxon][locus] =  msb_alignment.sequences
        self.assertEqual(alignments, self.expected_recoded_alignments)

        concat_alignments = {}
        for comparison_label, data in parsing.Loci.iter_from_msbayes_config(
            cfg,
            remove_triallelic_sites = False,
            convert_to_binary = False,
        ):
            s = StringIO()
            data.write_phylip(s)
            s.seek(0)
            sequences = []
            locus_2_idx = 0
            for i, line in enumerate(s):
                if i < 1:
                    continue
                line = line.strip()
                if line == '':
                    continue
                data = line.split()
                if len(data) == 2:
                    sequences.append(data)
                elif len(data) == 1:
                    sequences[locus_2_idx][1] += data[0]
                    locus_2_idx += 1
                else:
                    raise Exception('Unexpected phylip file line in test')
            concat_alignments[comparison_label] = sequences
        self.assertEqual(self.expected_concat_alignments_recode, concat_alignments)

    def test_parsing_dpp_msbayes_config(self):
        cfg = parsing.MsbayesConfig(
            self.dpp_cfg_path,
            recode_ambig_states_as_missing = False,
        )
        self.assertEqual(cfg.npairs, 2)

        expected_taxa = ["species1", "species2"]
        taxa = sorted(cfg.taxa)
        self.assertEqual(expected_taxa, taxa)

        self.assertEqual(cfg.using_dpp_model(), True)
        self.assertEqual(cfg.using_time_in_subs_per_site(), False)

        expected_eco_model = (
            'event_model_prior:\n'
            '    dirichlet_process:\n'
            '        parameters:\n'
            '            concentration:\n'
            '                estimate: true\n'
            '                prior:\n'
            '                    gamma_distribution:\n'
            '                        shape: 1.5\n'
            '                        scale: 5.0\n'
            '\n'
            'event_time_prior:\n'
            '    gamma_distribution:\n'
            '        shape: 1.0\n'
            '        scale: 0.05\n'
            '\n'
            'global_comparison_settings:\n'
            '    ploidy: 1\n'
            '    genotypes_are_diploid: false\n'
            '    markers_are_dominant: false\n'
            '    constant_sites_removed: false\n'
            '    population_name_delimiter: \" \"\n'
            '    population_name_is_prefix: false\n'
            '    equal_population_sizes: false\n'
            '    parameters:\n'
            '        freq_1:\n'
            '            value: 0.5\n'
            '            estimate: false\n'
            '        mutation_rate:\n'
            '            value: 1.0\n'
            '            estimate: false\n'
            '        root_relative_population_size:\n'
            '            value: 1.0\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 100.0\n'
            '                    scale: 0.01\n'
            '        population_size:\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 1.0\n'
            f'                    scale: {0.005/4.0}\n'
        )

        s = StringIO()
        cfg.write_ecoevolity_model_settings(s)
        self.assertEqual(expected_eco_model, s.getvalue())

        a = cfg.sample_table.alignments
        alignments = {}
        for taxon, locus_data in a.items():
            if not taxon in alignments:
                alignments[taxon] = {}
            for locus, msb_alignment in locus_data.items():
                alignments[taxon][locus] =  msb_alignment.sequences
        self.assertEqual(alignments, self.expected_alignments)

        concat_alignments = {}
        for comparison_label, data in parsing.Loci.iter_from_msbayes_config(
            cfg,
            remove_triallelic_sites = False,
            convert_to_binary = False,
        ):
            s = StringIO()
            data.write_phylip(s)
            s.seek(0)
            sequences = []
            locus_2_idx = 0
            for i, line in enumerate(s):
                if i < 1:
                    continue
                line = line.strip()
                if line == '':
                    continue
                data = line.split()
                if len(data) == 2:
                    sequences.append(data)
                elif len(data) == 1:
                    sequences[locus_2_idx][1] += data[0]
                    locus_2_idx += 1
                else:
                    raise Exception('Unexpected phylip file line in test')
            concat_alignments[comparison_label] = sequences
        self.assertEqual(self.expected_concat_alignments, concat_alignments)

    def test_parsing_dpp_msbayes_config_recode(self):
        cfg = parsing.MsbayesConfig(
            self.dpp_cfg_path,
            recode_ambig_states_as_missing = True,
        )
        self.assertEqual(cfg.npairs, 2)

        expected_taxa = ["species1", "species2"]
        taxa = sorted(cfg.taxa)
        self.assertEqual(expected_taxa, taxa)

        self.assertEqual(cfg.using_dpp_model(), True)
        self.assertEqual(cfg.using_time_in_subs_per_site(), False)

        expected_eco_model = (
            'event_model_prior:\n'
            '    dirichlet_process:\n'
            '        parameters:\n'
            '            concentration:\n'
            '                estimate: true\n'
            '                prior:\n'
            '                    gamma_distribution:\n'
            '                        shape: 1.5\n'
            '                        scale: 5.0\n'
            '\n'
            'event_time_prior:\n'
            '    gamma_distribution:\n'
            '        shape: 1.0\n'
            '        scale: 0.05\n'
            '\n'
            'global_comparison_settings:\n'
            '    ploidy: 1\n'
            '    genotypes_are_diploid: false\n'
            '    markers_are_dominant: false\n'
            '    constant_sites_removed: false\n'
            '    population_name_delimiter: \" \"\n'
            '    population_name_is_prefix: false\n'
            '    equal_population_sizes: false\n'
            '    parameters:\n'
            '        freq_1:\n'
            '            value: 0.5\n'
            '            estimate: false\n'
            '        mutation_rate:\n'
            '            value: 1.0\n'
            '            estimate: false\n'
            '        root_relative_population_size:\n'
            '            value: 1.0\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 100.0\n'
            '                    scale: 0.01\n'
            '        population_size:\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 1.0\n'
            f'                    scale: {0.005/4.0}\n'
        )

        s = StringIO()
        cfg.write_ecoevolity_model_settings(s)
        self.assertEqual(expected_eco_model, s.getvalue())

        a = cfg.sample_table.alignments
        alignments = {}
        for taxon, locus_data in a.items():
            if not taxon in alignments:
                alignments[taxon] = {}
            for locus, msb_alignment in locus_data.items():
                alignments[taxon][locus] =  msb_alignment.sequences
        self.assertEqual(alignments, self.expected_recoded_alignments)

        concat_alignments = {}
        for comparison_label, data in parsing.Loci.iter_from_msbayes_config(
            cfg,
            remove_triallelic_sites = False,
            convert_to_binary = False,
        ):
            s = StringIO()
            data.write_phylip(s)
            s.seek(0)
            sequences = []
            locus_2_idx = 0
            for i, line in enumerate(s):
                if i < 1:
                    continue
                line = line.strip()
                if line == '':
                    continue
                data = line.split()
                if len(data) == 2:
                    sequences.append(data)
                elif len(data) == 1:
                    sequences[locus_2_idx][1] += data[0]
                    locus_2_idx += 1
                else:
                    raise Exception('Unexpected phylip file line in test')
            concat_alignments[comparison_label] = sequences
        self.assertEqual(self.expected_concat_alignments_recode, concat_alignments)

    def test_parsing_dpp_sps_msbayes_config(self):
        cfg = parsing.MsbayesConfig(
            self.dpp_sps_cfg_path,
            recode_ambig_states_as_missing = False,
        )
        self.assertEqual(cfg.npairs, 2)

        expected_taxa = ["species1", "species2"]
        taxa = sorted(cfg.taxa)
        self.assertEqual(expected_taxa, taxa)

        self.assertEqual(cfg.using_dpp_model(), True)
        self.assertEqual(cfg.using_time_in_subs_per_site(), True)

        expected_eco_model = (
            'event_model_prior:\n'
            '    dirichlet_process:\n'
            '        parameters:\n'
            '            concentration:\n'
            '                estimate: true\n'
            '                prior:\n'
            '                    gamma_distribution:\n'
            '                        shape: 1.5\n'
            '                        scale: 5.0\n'
            '\n'
            'event_time_prior:\n'
            '    gamma_distribution:\n'
            '        shape: 1.0\n'
            '        scale: 0.1\n'
            '\n'
            'global_comparison_settings:\n'
            '    ploidy: 1\n'
            '    genotypes_are_diploid: false\n'
            '    markers_are_dominant: false\n'
            '    constant_sites_removed: false\n'
            '    population_name_delimiter: \" \"\n'
            '    population_name_is_prefix: false\n'
            '    equal_population_sizes: false\n'
            '    parameters:\n'
            '        freq_1:\n'
            '            value: 0.5\n'
            '            estimate: false\n'
            '        mutation_rate:\n'
            '            value: 1.0\n'
            '            estimate: false\n'
            '        root_relative_population_size:\n'
            '            value: 1.0\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 100.0\n'
            '                    scale: 0.01\n'
            '        population_size:\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 1.0\n'
            f'                    scale: {0.005/4.0}\n'
        )

        s = StringIO()
        cfg.write_ecoevolity_model_settings(s)
        self.assertEqual(expected_eco_model, s.getvalue())

        a = cfg.sample_table.alignments
        alignments = {}
        for taxon, locus_data in a.items():
            if not taxon in alignments:
                alignments[taxon] = {}
            for locus, msb_alignment in locus_data.items():
                alignments[taxon][locus] =  msb_alignment.sequences
        self.assertEqual(alignments, self.expected_alignments)

        concat_alignments = {}
        for comparison_label, data in parsing.Loci.iter_from_msbayes_config(
            cfg,
            remove_triallelic_sites = False,
            convert_to_binary = False,
        ):
            s = StringIO()
            data.write_phylip(s)
            s.seek(0)
            sequences = []
            locus_2_idx = 0
            for i, line in enumerate(s):
                if i < 1:
                    continue
                line = line.strip()
                if line == '':
                    continue
                data = line.split()
                if len(data) == 2:
                    sequences.append(data)
                elif len(data) == 1:
                    sequences[locus_2_idx][1] += data[0]
                    locus_2_idx += 1
                else:
                    raise Exception('Unexpected phylip file line in test')
            concat_alignments[comparison_label] = sequences
        self.assertEqual(self.expected_concat_alignments, concat_alignments)

    def test_parsing_dpp_sps_msbayes_config_recode(self):
        cfg = parsing.MsbayesConfig(
            self.dpp_sps_cfg_path,
            recode_ambig_states_as_missing = True,
        )
        self.assertEqual(cfg.npairs, 2)

        expected_taxa = ["species1", "species2"]
        taxa = sorted(cfg.taxa)
        self.assertEqual(expected_taxa, taxa)

        self.assertEqual(cfg.using_dpp_model(), True)
        self.assertEqual(cfg.using_time_in_subs_per_site(), True)

        expected_eco_model = (
            'event_model_prior:\n'
            '    dirichlet_process:\n'
            '        parameters:\n'
            '            concentration:\n'
            '                estimate: true\n'
            '                prior:\n'
            '                    gamma_distribution:\n'
            '                        shape: 1.5\n'
            '                        scale: 5.0\n'
            '\n'
            'event_time_prior:\n'
            '    gamma_distribution:\n'
            '        shape: 1.0\n'
            '        scale: 0.1\n'
            '\n'
            'global_comparison_settings:\n'
            '    ploidy: 1\n'
            '    genotypes_are_diploid: false\n'
            '    markers_are_dominant: false\n'
            '    constant_sites_removed: false\n'
            '    population_name_delimiter: \" \"\n'
            '    population_name_is_prefix: false\n'
            '    equal_population_sizes: false\n'
            '    parameters:\n'
            '        freq_1:\n'
            '            value: 0.5\n'
            '            estimate: false\n'
            '        mutation_rate:\n'
            '            value: 1.0\n'
            '            estimate: false\n'
            '        root_relative_population_size:\n'
            '            value: 1.0\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 100.0\n'
            '                    scale: 0.01\n'
            '        population_size:\n'
            '            estimate: true\n'
            '            prior:\n'
            '                gamma_distribution:\n'
            '                    shape: 1.0\n'
            f'                    scale: {0.005/4.0}\n'
        )

        s = StringIO()
        cfg.write_ecoevolity_model_settings(s)
        self.assertEqual(expected_eco_model, s.getvalue())

        a = cfg.sample_table.alignments
        alignments = {}
        for taxon, locus_data in a.items():
            if not taxon in alignments:
                alignments[taxon] = {}
            for locus, msb_alignment in locus_data.items():
                alignments[taxon][locus] =  msb_alignment.sequences
        self.assertEqual(alignments, self.expected_recoded_alignments)

        concat_alignments = {}
        for comparison_label, data in parsing.Loci.iter_from_msbayes_config(
            cfg,
            remove_triallelic_sites = False,
            convert_to_binary = False,
        ):
            s = StringIO()
            data.write_phylip(s)
            s.seek(0)
            sequences = []
            locus_2_idx = 0
            for i, line in enumerate(s):
                if i < 1:
                    continue
                line = line.strip()
                if line == '':
                    continue
                data = line.split()
                if len(data) == 2:
                    sequences.append(data)
                elif len(data) == 1:
                    sequences[locus_2_idx][1] += data[0]
                    locus_2_idx += 1
                else:
                    raise Exception('Unexpected phylip file line in test')
            concat_alignments[comparison_label] = sequences
        self.assertEqual(self.expected_concat_alignments_recode, concat_alignments)


if __name__ == '__main__':
    unittest.main()

