#! /usr/bin/env python

import unittest
import os
import sys
import logging
import io
from contextlib import redirect_stdout

from pycoevolity.cli import loci2alignment
from pycoevolity.test import TestLevel
from pycoevolity.test.support.pycoevolity_test_case import PycoevolityTestCase

_LOG = logging.getLogger(__name__)


class Loci2AlignmentTestCase(PycoevolityTestCase):
    def setUp(self):
        self.set_up()

    def tearDown(self):
        self.tear_down()

    def test_2_loci(self):
        loci_path = self.get_test_path(prefix = "temp-loci")
        loci_str = """CDS_4485_Cyrtodactylus_annulatus_Bohol           TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
CDS_4658_Cyrtodactylus_annulatus_Bohol           TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
CDS_4659_Cyrtodactylus_annulatus_Bohol           TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
CDS_4660_Cyrtodactylus_annulatus_Bohol           TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
RMB_8043_Cyrtodactylus_annulatus_CamiguinSur     TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTMTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
RMB_8201_Cyrtodactylus_annulatus_CamiguinSur     TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
RMB_8220_Cyrtodactylus_annulatus_CamiguinSur     TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
RMB_8232_Cyrtodactylus_annulatus_CamiguinSur     TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTMTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
//                                                                                                        *                                 |0|
CDS_4485_Cyrtodactylus_annulatus_Bohol           GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
CDS_4659_Cyrtodactylus_annulatus_Bohol           --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAANTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
CDS_4660_Cyrtodactylus_annulatus_Bohol           --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTNTTTCTGCAAAACNATTTTCAGCATGCTCTCTACCATCAACCTTCACA
RMB_8043_Cyrtodactylus_annulatus_CamiguinSur     GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
RMB_8201_Cyrtodactylus_annulatus_CamiguinSur     GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
RMB_8220_Cyrtodactylus_annulatus_CamiguinSur     GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
RMB_8232_Cyrtodactylus_annulatus_CamiguinSur     --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
//                                                                                                                                                      |1|
"""
        with open(loci_path, "w") as out:
            out.write(loci_str)
        args = [
            "--prefix", "foo-",
            "--suffix=-bar",
            "--charsets",
            loci_path,
        ]
        out_stream = io.StringIO()
        with redirect_stdout(out_stream):
            loci2alignment.main_nexus(argv = args)
        out_str = out_stream.getvalue()
        
        expected_out_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=8;
    TAXLABELS
        foo-CDS_4485_Cyrtodactylus_annulatus_Bohol-bar
        foo-CDS_4658_Cyrtodactylus_annulatus_Bohol-bar
        foo-CDS_4659_Cyrtodactylus_annulatus_Bohol-bar
        foo-CDS_4660_Cyrtodactylus_annulatus_Bohol-bar
        foo-RMB_8043_Cyrtodactylus_annulatus_CamiguinSur-bar
        foo-RMB_8201_Cyrtodactylus_annulatus_CamiguinSur-bar
        foo-RMB_8220_Cyrtodactylus_annulatus_CamiguinSur-bar
        foo-RMB_8232_Cyrtodactylus_annulatus_CamiguinSur-bar
    ;
END;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR=194;
    FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=YES;
    MATRIX
        foo-CDS_4485_Cyrtodactylus_annulatus_Bohol-bar          TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-CDS_4658_Cyrtodactylus_annulatus_Bohol-bar          TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-CDS_4659_Cyrtodactylus_annulatus_Bohol-bar          TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-CDS_4660_Cyrtodactylus_annulatus_Bohol-bar          TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-RMB_8043_Cyrtodactylus_annulatus_CamiguinSur-bar    TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTMTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-RMB_8201_Cyrtodactylus_annulatus_CamiguinSur-bar    TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-RMB_8220_Cyrtodactylus_annulatus_CamiguinSur-bar    TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-RMB_8232_Cyrtodactylus_annulatus_CamiguinSur-bar    TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTMTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT

        foo-CDS_4485_Cyrtodactylus_annulatus_Bohol-bar          GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-CDS_4658_Cyrtodactylus_annulatus_Bohol-bar          ???????????????????????????????????????????????????????????????????????????????????????????????????????
        foo-CDS_4659_Cyrtodactylus_annulatus_Bohol-bar          --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAA?TATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-CDS_4660_Cyrtodactylus_annulatus_Bohol-bar          --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTT?TTTCTGCAAAAC?ATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-RMB_8043_Cyrtodactylus_annulatus_CamiguinSur-bar    GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-RMB_8201_Cyrtodactylus_annulatus_CamiguinSur-bar    GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-RMB_8220_Cyrtodactylus_annulatus_CamiguinSur-bar    GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-RMB_8232_Cyrtodactylus_annulatus_CamiguinSur-bar    --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
    ;
END;

BEGIN SETS;
    CHARSET locus1=1-91;
    CHARSET locus2=92-194;
END;
"""
        # lines = out_str.split("\n")
        # elines = expected_out_str.split("\n")
        # self.assertEqual(len(lines), len(elines))
        # for i in range(len(lines)):
        #     if not lines[i] == elines[i]:
        #         print(lines[i])
        #         print(elines[i])
        #     self.assertEqual(lines[i], elines[i])
        self.assertEqual(out_str, expected_out_str)

    def test_2_loci_remove_sample(self):
        loci_path = self.get_test_path(prefix = "temp-loci")
        loci_str = """CDS_4485_Cyrtodactylus_annulatus_Bohol           TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
CDS_4658_Cyrtodactylus_annulatus_Bohol           TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
CDS_4659_Cyrtodactylus_annulatus_Bohol           TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
CDS_4660_Cyrtodactylus_annulatus_Bohol           TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
RMB_8043_Cyrtodactylus_annulatus_CamiguinSur     TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTMTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
RMB_8201_Cyrtodactylus_annulatus_CamiguinSur     TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
RMB_8220_Cyrtodactylus_annulatus_CamiguinSur     TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
RMB_8232_Cyrtodactylus_annulatus_CamiguinSur     TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTMTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
//                                                                                                        *                                 |0|
CDS_4485_Cyrtodactylus_annulatus_Bohol           GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
CDS_4659_Cyrtodactylus_annulatus_Bohol           --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAANTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
CDS_4660_Cyrtodactylus_annulatus_Bohol           --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTNTTTCTGCAAAACNATTTTCAGCATGCTCTCTACCATCAACCTTCACA
RMB_8043_Cyrtodactylus_annulatus_CamiguinSur     GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
RMB_8201_Cyrtodactylus_annulatus_CamiguinSur     GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
RMB_8220_Cyrtodactylus_annulatus_CamiguinSur     GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
RMB_8232_Cyrtodactylus_annulatus_CamiguinSur     --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
//                                                                                                                                                      |1|
"""
        with open(loci_path, "w") as out:
            out.write(loci_str)
        args = [
            "--prefix", "foo-",
            "--suffix=-bar",
            "--charsets",
            "--sample-to-delete", "CDS_4658_Cyrtodactylus_annulatus_Bohol",
            loci_path,
        ]
        out_stream = io.StringIO()
        with redirect_stdout(out_stream):
            loci2alignment.main_nexus(argv = args)
        out_str = out_stream.getvalue()
        
        expected_out_str = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=7;
    TAXLABELS
        foo-CDS_4485_Cyrtodactylus_annulatus_Bohol-bar
        foo-CDS_4659_Cyrtodactylus_annulatus_Bohol-bar
        foo-CDS_4660_Cyrtodactylus_annulatus_Bohol-bar
        foo-RMB_8043_Cyrtodactylus_annulatus_CamiguinSur-bar
        foo-RMB_8201_Cyrtodactylus_annulatus_CamiguinSur-bar
        foo-RMB_8220_Cyrtodactylus_annulatus_CamiguinSur-bar
        foo-RMB_8232_Cyrtodactylus_annulatus_CamiguinSur-bar
    ;
END;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR=194;
    FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=YES;
    MATRIX
        foo-CDS_4485_Cyrtodactylus_annulatus_Bohol-bar          TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-CDS_4659_Cyrtodactylus_annulatus_Bohol-bar          TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-CDS_4660_Cyrtodactylus_annulatus_Bohol-bar          TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-RMB_8043_Cyrtodactylus_annulatus_CamiguinSur-bar    TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTMTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-RMB_8201_Cyrtodactylus_annulatus_CamiguinSur-bar    TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-RMB_8220_Cyrtodactylus_annulatus_CamiguinSur-bar    TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTCTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT
        foo-RMB_8232_Cyrtodactylus_annulatus_CamiguinSur-bar    TCTTTTGACAAAGCTGCAAAGCACAATCTTCCGTATAAGACTCACATGATGCCCATTMTTTCAGTAAGATTTGTGCCAGAAAAGTTTTCAT

        foo-CDS_4485_Cyrtodactylus_annulatus_Bohol-bar          GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-CDS_4659_Cyrtodactylus_annulatus_Bohol-bar          --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAA?TATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-CDS_4660_Cyrtodactylus_annulatus_Bohol-bar          --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTT?TTTCTGCAAAAC?ATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-RMB_8043_Cyrtodactylus_annulatus_CamiguinSur-bar    GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-RMB_8201_Cyrtodactylus_annulatus_CamiguinSur-bar    GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-RMB_8220_Cyrtodactylus_annulatus_CamiguinSur-bar    GTATTGCATATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
        foo-RMB_8232_Cyrtodactylus_annulatus_CamiguinSur-bar    --------TATGCACTATAATCTAGTAAATTTCAGAGTGAAATTATCAGGACCTTCTTTCTGCAAAACTATTTTCAGCATGCTCTCTACCATCAACCTTCACA
    ;
END;

BEGIN SETS;
    CHARSET locus1=1-91;
    CHARSET locus2=92-194;
END;
"""
        # lines = out_str.split("\n")
        # elines = expected_out_str.split("\n")
        # self.assertEqual(len(lines), len(elines))
        # for i in range(len(lines)):
        #     if not lines[i] == elines[i]:
        #         print(lines[i])
        #         print(elines[i])
        #     self.assertEqual(lines[i], elines[i])
        self.assertEqual(out_str, expected_out_str)

if __name__ == '__main__':
    unittest.main()

