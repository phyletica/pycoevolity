#! /usr/bin/env python

"""
CLI program for generating dummy biallelic data files.
"""

import os
import sys
import argparse

import pycoevolity


def main(argv = sys.argv):
    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-s', '--nspecies',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_int,
            default = 2,
            help = ('The number of species or populations.'))
    parser.add_argument('-g', '--ngenomes',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_int,
            default = 10,
            help = ('The number of genomes sampled per population.'))
    parser.add_argument('-c', '--ncharacters',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_int,
            default = 1000,
            help = ('The number of biallelic characters.'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            default = "sp",
            help = ('Prefix for species labels.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    pycoevolity.write_dummy_biallelic_data_file(
            nspecies = args.nspecies,
            ngenomes = args.ngenomes,
            ncharacters = args.ncharacters,
            prefix = args.prefix,
            out = sys.stdout)


if __name__ == "__main__":
    main()
