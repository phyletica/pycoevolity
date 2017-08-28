#! /usr/bin/env python

"""
CLI program for converting a *.loci file from ipyrad to nexus format.
"""

import os
import sys
import argparse

import sumcoevolity

def main(argv = sys.argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('loci_path',
            metavar = 'IPYRAD-LOCI-FILE-PATH',
            type = sumcoevolity.argparse_utils.arg_is_file,
            help = ('Path to the ipyrad loci file.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    data = sumcoevolity.parsing.PyradLoci(args.loci_path)
    sys.stderr.write("Number of taxa:  {0}\n".format(data.number_of_taxa))
    sys.stderr.write("Number of loci:  {0}\n".format(data.number_of_loci))
    sys.stderr.write("Number of sites: {0}\n".format(data.number_of_sites))
    data.write_nexus()

if __name__ == "__main__":
    main()
