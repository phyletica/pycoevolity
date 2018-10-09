#! /usr/bin/env python

"""
CLI program for converting a *.loci file from ipyrad to dpp-msbayes format.
"""

import os
import sys
import random
import argparse

import pycoevolity

def get_dpp_msbayes_preamble():
    return """concentrationShape = 10000.0
concentrationScale = 0.0001414216
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

def main(argv = sys.argv):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('loci_paths',
            metavar = 'IPYRAD-LOCI-FILE-PATH',
            nargs = '+',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Paths to the ipyrad loci file.'))
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

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    rng = random.Random()
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    rng.seed(args.seed)

    sys.stderr.write("Command: {0}\n\n".format(" ".join(argv)))

    path_suffix = ""
    if args.number_of_loci:
        path_suffix = "-{0}-{1}".format(args.number_of_loci, args.seed)

    if not args.output_dir:
        args.output_dir = os.path.join(
                os.getcwd(),
                "loci2dppmsbayes-output{0}".format(path_suffix)
                )
        pycoevolity.fileio.make_directory(args.output_dir)

    config_path = os.path.join(
            args.output_dir,
            "dpp-msbayes{0}.cfg".format(path_suffix)
            )
    if os.path.exists(config_path):
        raise Exception("The path {0!r} already exists. "
                        "Please designate a different "
                        "directory.".format(config_path))

    fasta_dir = os.path.join(args.output_dir, "fasta-alignments")
    nexus_dir = os.path.join(args.output_dir, "nexus-alignments")
    make_directory(fasta_dir)
    make_directory(nexus_dir)
    relative_fasta_dir = os.path.relpath(fasta_dir, start = args.output_dir)

    pair_number_buffer = len(str(len(args.loci_paths)))

    with open(config_path, 'w') as out_stream:
        out_stream.write(get_dpp_msbayes_preamble())
        out_stream.write("\nBEGIN SAMPLE_TBL\n")

        for i, loci_path in enumerate(args.loci_paths):
            sys.stderr.write("Parsing {0!r}...\n".format(loci_path))
            data = pycoevolity.parsing.PyradLoci(loci_path,
                    remove_triallelic_sites = args.remove_triallelic_sites,
                    sequence_ids_to_remove = args.sample_to_delete)
            pair_label = "{pair_num:0{buffer_size}d}".format(
                    pair_num = i,
                    buffer_size = pair_number_buffer)
            if args.population_name_is_prefix:
                data.label_prefix = pair_label
            else:
                data.label_suffix = pair_label
            seqs_removed = data.removed_sequence_counts
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

            data.write_fasta_files(
                    directory = fasta_dir,
                    pair_label = pair_label,
                    population_name_delimiter = args.population_name_delimiter,
                    population_name_is_prefix = args.population_name_is_prefix,
                    write_sample_table = True,
                    sample_table_stream = out_stream,
                    sample_table_directory = os.path.dirname(config_path))

            nexus_prefix = os.path.splitext(os.path.basename(loci_path))[0]
            if loci_path.endswith("gz"):
                nexus_prefix = os.path.splitext(os.path.basename(nexus_prefix))[0]
            nexus_path = os.path.join(
                    nexus_dir,
                    "{0}{1}.nex".format(
                            nexus_prefix,
                            path_suffix)
                    )
            if os.path.exists(nexus_path):
                raise Exception("The path {0!r} already exists. "
                                "Please designate a different "
                                "directory.".format(nexus_path))
            with open(nexus_path, 'w') as nexus_stream:
                data.write_nexus(stream = nexus_stream)
        out_stream.write("\nEND SAMPLE_TBL\n")

if __name__ == "__main__":
    main()
