#! /usr/bin/env python

"""
CLI program for converting a *.loci file from ipyrad to nexus format.
"""

import os
import sys
import argparse

import pycoevolity

def main(argv = sys.argv, write_method = "write_nexus"):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('loci_path',
            metavar = 'IPYRAD-LOCI-FILE-PATH',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Path to the ipyrad loci file.'))
    parser.add_argument('-d', '--sample-to-delete',
            action = 'append',
            type = str,
            metavar = "SAMPLE-ID",
            help = ('The ID of a sequence that should not be included in the '
                    'nexus alignment. This option can be used multiple times '
                    'to remove multiple samples.'))
    parser.add_argument('-r', '--remove-triallelic-sites',
            action = 'store_true',
            help = ('By default, sites with more than two alleles are left in '
                    'the alignment, as is. When this option is specified, '
                    'these sites are removed.'))
    parser.add_argument('-b', '--convert-to-binary',
            action = 'store_true',
            help = ('By default, DNA characters are maintained in the output '
                    'alignment. When this option is specified, the characters '
                    'are converted to binary states. Ecoevolity will accept '
                    'DNA characters and convert to binary internally; this '
                    'option is primarily for testing.'))
    parser.add_argument('-p', '--prefix',
            type = str,
            metavar = "SEQUENCE-LABEL-PREFIX",
            help = ('A prefix to append to every sequence label. This can be '
                    'useful for ensuring that all population labels are '
                    'unique across pairs.'))
    parser.add_argument('-s', '--suffix',
            type = str,
            metavar = "SEQUENCE-LABEL-SUFFIX",
            help = ('A suffix to append to every sequence label. This can be '
                    'useful for ensuring that all population labels are '
                    'unique across pairs.'))
    parser.add_argument('-c', '--charsets',
            action = 'store_true',
            help = ('Include charsets block in output nexus file. This option '
                    'is ignored if output format is not nexus.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    data = pycoevolity.parsing.PyradLoci(args.loci_path,
            remove_triallelic_sites = args.remove_triallelic_sites,
            convert_to_binary = args.convert_to_binary,
            sequence_ids_to_remove = args.sample_to_delete)
    if args.prefix:
        data.label_prefix = args.prefix
    if args.suffix:
        data.label_suffix = args.suffix
    seqs_removed = data.removed_sequence_counts
    sys.stderr.write("Command: {0}\n".format(" ".join(argv)))
    sys.stderr.write("\tNumber of taxa:  {0}\n".format(data.number_of_taxa))
    sys.stderr.write("\tNumber of loci:  {0}\n".format(data.number_of_loci))
    sys.stderr.write("\tNumber of sites: {0}\n".format(data.number_of_sites))
    sys.stderr.write("\tNumber of triallelic sites found:   {0}\n".format(
            data.number_of_triallelic_sites_found))
    sys.stderr.write("\tNumber of triallelic sites removed: {0}\n".format(
            data.number_of_triallelic_sites_removed))
    if seqs_removed:
        sys.stderr.write("\tNumber of sequences removed:\n")
        for n, c in seqs_removed.items():
            sys.stderr.write("\t\t{0}: {1}\n".format(n, c))
            
    write_kwargs = {}
    if write_method == "write_nexus":
        if args.charsets:
            write_kwargs["include_charset_block"] = True
    getattr(data, write_method)(**write_kwargs)

def main_nexus(argv = sys.argv):
    main(argv = argv, write_method = "write_nexus")

def main_phylip(argv = sys.argv):
    main(argv = argv, write_method = "write_phylip")

if __name__ == "__main__":
    main()
