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

if __name__ == '__main__':
    unittest.main()

