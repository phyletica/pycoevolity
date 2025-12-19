#! /usr/bin/env python

import unittest
import os
import sys
import logging
import io
from contextlib import redirect_stdout

from pycoevolity.cli import fastas2alignment
from pycoevolity.test import TestLevel
from pycoevolity.test.support.pycoevolity_test_case import PycoevolityTestCase

_LOG = logging.getLogger(__name__)


class Fastas2AlignmentTestCase(PycoevolityTestCase):
    def setUp(self):
        self.set_up()

    def tearDown(self):
        self.tear_down()

    def test_2_loci(self):
        tmp_fasta_1 = self.get_test_path(prefix = "temp-fasta-1")
        tmp_fasta_2 = self.get_test_path(prefix = "temp-fasta-2")
        fasta_1 = """
>seq-a
acgt-n?
>seq-b
acgc-at
"""
        fasta_2 = """
>seq-b
GGTC-A
>seq-c
GGCTTT
"""
        with open(tmp_fasta_1, "w") as out:
            out.write(fasta_1)
        with open(tmp_fasta_2, "w") as out:
            out.write(fasta_2)
        args = [
            "--treat-n-as-missing",
            "--prefix", "foo-",
            "--suffix=-bar",
            "--charsets",
            tmp_fasta_1,
            tmp_fasta_2,
        ]
        out_stream = io.StringIO()
        with redirect_stdout(out_stream):
            fastas2alignment.main_nexus(argv = args)
        out_str = out_stream.getvalue()
        
        expected_out_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=3;
    TAXLABELS
        foo-seq-a-bar
        foo-seq-b-bar
        foo-seq-c-bar
    ;
END;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR=13;
    FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=YES;
    MATRIX
        foo-seq-a-bar    ACGT-??
        foo-seq-b-bar    ACGC-AT
        foo-seq-c-bar    ???????

        foo-seq-a-bar    ??????
        foo-seq-b-bar    GGTC-A
        foo-seq-c-bar    GGCTTT
    ;
END;

BEGIN SETS;
    CHARSET locus1=1-7;
    CHARSET locus2=8-13;
END;
"""
        self.assertEqual(out_str, expected_out_str)

    def test_ambig_char_error(self):
        tmp_fasta_1 = self.get_test_path(prefix = "temp-fasta-1")
        tmp_fasta_2 = self.get_test_path(prefix = "temp-fasta-2")
        fasta_1 = """
>seq-a
acgt-n?
>seq-b
acgc-at
"""
        fasta_2 = """
>seq-b
GGTC-A
>seq-c
GGCTTT
"""
        with open(tmp_fasta_1, "w") as out:
            out.write(fasta_1)
        with open(tmp_fasta_2, "w") as out:
            out.write(fasta_2)
        args = [
            "--prefix", "foo-",
            "--suffix=-bar",
            "--charsets",
            tmp_fasta_1,
            tmp_fasta_2,
        ]
        out_stream = io.StringIO()
        with redirect_stdout(out_stream):
            self.assertRaises(Exception, fastas2alignment.main_nexus, *args)
        out_str = out_stream.getvalue()

    def test_2_loci_remove_sample(self):
        tmp_fasta_1 = self.get_test_path(prefix = "temp-fasta-1")
        tmp_fasta_2 = self.get_test_path(prefix = "temp-fasta-2")
        fasta_1 = """
>seq-a
acgt-n?
>seq-b
acgc-at
"""
        fasta_2 = """
>seq-b
GGTC-A
>seq-c
GGCTTT
"""
        with open(tmp_fasta_1, "w") as out:
            out.write(fasta_1)
        with open(tmp_fasta_2, "w") as out:
            out.write(fasta_2)
        args = [
            "--treat-n-as-missing",
            "--prefix", "foo-",
            "--suffix=-bar",
            "--charsets",
            "--sample-to-delete", "seq-b",
            tmp_fasta_1,
            tmp_fasta_2,
        ]
        out_stream = io.StringIO()
        with redirect_stdout(out_stream):
            fastas2alignment.main_nexus(argv = args)
        out_str = out_stream.getvalue()
        
        expected_out_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=2;
    TAXLABELS
        foo-seq-a-bar
        foo-seq-c-bar
    ;
END;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR=13;
    FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=YES;
    MATRIX
        foo-seq-a-bar    ACGT-??
        foo-seq-c-bar    ???????

        foo-seq-a-bar    ??????
        foo-seq-c-bar    GGCTTT
    ;
END;

BEGIN SETS;
    CHARSET locus1=1-7;
    CHARSET locus2=8-13;
END;
"""
        self.assertEqual(out_str, expected_out_str)
        

if __name__ == '__main__':
    unittest.main()

