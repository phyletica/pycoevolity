#! /usr/bin/env python

"""
CLI program for converting fasta files to fasta files that all contain the
union of all sequence labels.
"""

import os
import sys
import random
import argparse

import pycoevolity

def main(argv = sys.argv):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('fasta_paths',
            metavar = 'FASTA-FILE-PATH',
            nargs = '+',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Paths to fasta files.'))
    parser.add_argument('-o', '--output-dir',
            metavar = 'OUTPUT-DIRECTORY',
            type = pycoevolity.argparse_utils.arg_is_dir,
            help = ('Path to directory in which to write output files. '
                    'Default: \'pyco-fasta-files\' directory within the '
                    'current working directory.'))
    parser.add_argument('-n', '--number-of-loci',
            type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
            help = ('The number of loci to randomly sample. Default behavior '
                    'is to convert and output all the loci. If this option is '
                    'used, loci are randomly sampled without replacement, so '
                    'the number can not be larger than the number of loci. You '
                    'can use in conjuction with the \'--with-replacement\' '
                    'option if you want to randomly sample loci with '
                    'replacement.'))
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
    parser.add_argument('--seed',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_int,
            help = ('Seed for random number generator.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    rng = random.Random()
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    rng.seed(args.seed)

    if not args.output_dir:
        args.output_dir = os.path.join(
                os.getcwd(),
                "pyco-fasta-files"
                )
        pycoevolity.fileio.make_directory(args.output_dir)

    data = pycoevolity.parsing.Loci.from_fastas(args.fasta_paths,
            remove_triallelic_sites = args.remove_triallelic_sites,
            convert_to_binary = args.convert_to_binary,
            sequence_ids_to_remove = args.sample_to_delete,
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

    if args.number_of_loci:
        data.sample_loci(
                rng = rng,
                number_of_samples = args.number_of_loci,
                with_replacement = args.with_replacement)

    data.write_union_fasta_files(
            directory = args.output_dir)

if __name__ == "__main__":
    main()
