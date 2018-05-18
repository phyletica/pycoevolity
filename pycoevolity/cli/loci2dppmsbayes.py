#! /usr/bin/env python

"""
CLI program for converting a *.loci file from ipyrad to dpp-msbayes format.
"""

import os
import sys
import errno
import random
import argparse

import pycoevolity

def get_dpp_msbayes_preamble():
    return """concentrationShape = 10000.0
concentrationScale = 0.000182147
thetaShape = 5.0
thetaScale = 0.0004
ancestralThetaShape = 0
ancestralThetaScale = 0
thetaParameters = 012
tauShape = 1.0
tauScale = 0.1
timeInSubsPerSite = 1
bottleProportionShapeA = 0.0
bottleProportionShapeB = 0.0
bottleProportionShared = 0
numTauClasses = 0
"""

def make_directory(path):
    """
    Creates directory `path`, but suppresses error if `path` already exists.
    """
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise e

def main(argv = sys.argv):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('loci_path',
            metavar = 'IPYRAD-LOCI-FILE-PATH',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Path to the ipyrad loci file.'))
    parser.add_argument('--output-dir',
            type = pycoevolity.argparse_utils.arg_is_dir,
            help = ('Path to a directory where you want all the output to go. '
                    'Default: \'loci2dppmsbayes-output\' directory within the '
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
    parser.add_argument('--with-replacement',
            action = 'store_true',
            help = ('If \'--number-of-loci\' is specified, that number of '
                    'loci will be randomly subsampled without replacement. '
                    'This option will change that behavior to randomly sample '
                    'with replacement. If the \'--number-of-loci\' option is '
                    'NOT specified, this argument is ignored, and all loci '
                    'are converted and output.'))
    parser.add_argument('--population-name-delimiter',
            type = str,
            default = '_',
            help = ('The character used to delimit the population name from '
                    'the rest of the sequence label. '
                    'Default: \'_\''))
    parser.add_argument('--population-name-is-prefix',
            action = 'store_true',
            help = ('By default, the program looks for the population name '
                    'at the end of each sequence label, delimited by the '
                    '\'--population-name-delimiter]\'. This option tells  the '
                    'program to look at the beginning of each sequence label.'))
    parser.add_argument('--seed',
            action = 'store',
            type = pycoevolity.argparse_utils.arg_is_positive_int,
            help = ('Seed for random number generator.'))
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

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    if not args.output_dir:
        args.output_dir = os.path.join(os.getcwd(), "loci2dppmsbayes-output")
        make_directory(args.output_dir)

    rng = random.Random()
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    rng.seed(args.seed)

    data = pycoevolity.parsing.PyradLoci(args.loci_path,
            remove_triallelic_sites = args.remove_triallelic_sites,
            sequence_ids_to_remove = args.sample_to_delete)
    if args.prefix:
        data.label_prefix = args.prefix
    if args.suffix:
        data.label_suffix = args.suffix
    seqs_removed = data.removed_sequence_counts
    sys.stderr.write("Command: {0}\n".format(" ".join(argv)))
    sys.stderr.write("Stats on total data parsed:\n".format(" ".join(argv)))
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

    config_path = os.path.join(args.output_dir, "dpp-msbayes.cfg")
    if os.path.exists(config_path):
        raise Exception("The path {0!r} already exists. "
                        "Please designate a different directory.".format(
                                config_path))

    with open(config_path, 'w') as out_stream:
        out_stream.write(get_dpp_msbayes_preamble())
        out_stream.write("\nBEGIN SAMPLE_TBL\n")
        data.write_fasta_files(
                directory = args.output_dir,
                write_sample_table = True,
                population_name_delimiter = args.population_name_delimiter,
                population_name_is_prefix = args.population_name_is_prefix,
                stream = out_stream)
        out_stream.write("\nEND SAMPLE_TBL\n")
            
    data.write_nexus(stream = sys.stdout)

if __name__ == "__main__":
    main()
