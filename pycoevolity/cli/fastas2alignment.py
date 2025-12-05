#! /usr/bin/env python

"""
CLI program for converting fasta files to concatenated nexus format.
"""

import os
import sys
import argparse

import pycoevolity

def main(argv = sys.argv, write_method = "write_nexus"):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('fasta_paths',
            metavar = 'FASTA-FILE-PATH',
            nargs = '+',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Paths to fasta files.'))
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
    parser.add_argument('--treat-n-as-missing',
            action = 'store_true',
            help = (
                'By default, ambiguity codes are treated as polymorphisms, '
                'and symbols that encode more than two states cause an error, '
                'because they are not biallelic. When this option is specified '
                '\'N\'s will be treated as missing data (\'?\') and will not '
                'trigger this error.'
            ))
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
    parser.add_argument('-m', '--label-change-map',
            type = pycoevolity.argparse_utils.arg_is_file,
            metavar = "SEQ-LABEL-CHANGE-CSV-PATH",
            help = ('Path to a CSV file with sequence labels to find and '
                    'replace. The CSV file should have 2 columns: current '
                    'sequence labels in the first, and the replacements in '
                    'the second.'))
    parser.add_argument('-c', '--charsets',
            action = 'store_true',
            help = ('Include charsets block in output nexus file. This option '
                    'is ignored if output format is not nexus.'))
    parser.add_argument('--split',
            action = 'store_true',
            help = ('Randomly split loci into two output alignments. This '
                    'option should only be used for testing whether ecoevolity '
                    'will estimate co-divergence between the random sets of '
                    'loci from the same taxon.'))
    parser.add_argument('--subsample',
            type = int,
            default = 0,
            metavar = "NUMBER-OF-LOCI",
            help = ('Randomly subsample this number of loci without replacement.'))
    parser.add_argument('--seed',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_int,
            help = ('Seed for random number generator. This is only used for '
                    'the \'--split\' or \'--subsample\'  options.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    if args.split and (args.subsample > 0):
        msg = "ERROR: '--split' and '--subsample' cannot be used together"
        raise Exception(msg)

    rng = random.Random()
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    rng.seed(args.seed)

    label_change_map_path = None
    if args.label_change_map:
        label_change_map_path = args.label_change_map

    data = pycoevolity.parsing.Loci.from_fastas(args.fasta_paths,
            remove_triallelic_sites = args.remove_triallelic_sites,
            convert_to_binary = args.convert_to_binary,
            sequence_ids_to_remove = args.sample_to_delete,
            label_change_map_path = label_change_map_path,
            treat_n_as_missing = args.treat_n_as_missing)
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

    if args.split:
        data2 = data.split_loci(rng = rng,
                auto_annotate_labels = True)
        sys.stderr.write("\tNumber of loci in set 1: {0}\n".format(data.number_of_loci))
        sys.stderr.write("\tNumber of loci in set 2: {0}\n".format(data2.number_of_loci))
        getattr(data2, write_method)(**write_kwargs)
    elif args.subsample > 0:
        data.sample_loci(rng = rng,
                number_of_samples = args.subsample,
                with_replacement = False)
    getattr(data, write_method)(**write_kwargs)

def main_nexus(argv = sys.argv):
    main(argv = argv, write_method = "write_nexus")

def main_phylip(argv = sys.argv):
    main(argv = argv, write_method = "write_phylip")

if __name__ == "__main__":
    main()
