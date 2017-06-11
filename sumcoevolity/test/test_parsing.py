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

from sumcoevolity import parsing
from sumcoevolity.test.support.sumcoevolity_test_case import SumcoevolityTestCase
from sumcoevolity.test.support import package_paths

_LOG = logging.getLogger(__name__)

class GetDictFromSpreadsheetTestCase(SumcoevolityTestCase):

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

class EcoevolityStdOutTestCase(SumcoevolityTestCase):
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


if __name__ == '__main__':
    unittest.main()

