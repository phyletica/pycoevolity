#! /usr/bin/env python

"""
CLI program for converting MsBayes-like fasta files to concatenated nexus
format.
"""

import os
import sys
import argparse

import pycoevolity

def main(argv = sys.argv, write_method = "write_nexus"):
    pycoevolity.write_splash(sys.stderr)
    parser = argparse.ArgumentParser()

    parser.add_argument('msbayes_config_path',
            metavar = 'MSBAYES-CONFIG-PATH',
            type = pycoevolity.argparse_utils.arg_is_file,
            help = ('Path to MsBayes-like config file.'))
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
    parser.add_argument('-m', '--recode-ambig-states-as-missing',
            action = 'store_true',
            help = ('By default, ambiguous DNA characters are maintained in '
                    'the output alignment. When this option is specified, the '
                    'ambiguous characters are converted to missing (\'?\'). '
                    'This is useful when the sequences are haplotypes '
                    '(i.e., genotypes are haploid) and thus ambiguous '
                    'states do NOT represent heterozygous alleles.'))
    parser.add_argument('-c', '--charsets',
            action = 'store_true',
            help = ('Include charsets block in output nexus file. This option '
                    'is ignored if output format is not nexus.'))
    parser.add_argument('-p', '--comparison-prefix',
            metavar = 'COMPARISON-NEXUS-FILE-PREFIX',
            type = str,
            default = '',
            help = ('In the ecoevolity yaml config written to stdout, this '
                    'prefix will be added to the nexus file name for each '
                    'comparison.'))
    parser.add_argument('-o', '--output-dir',
            metavar = 'OUTPUT-DIRECTORY',
            type = pycoevolity.argparse_utils.arg_is_dir,
            help = ('Path to directory in which to write output files. '
                    'Default: \'pyco-nexus-files\' directory within the '
                    'same directory as the input config file.'))
    parser.add_argument('-y', '--yaml-config-only',
            action = 'store_true',
            help = ('Only write yaml config to stdout (i.e., do not produce '
                    'nexus files.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    if not args.output_dir:
        args.output_dir = os.path.join(
            os.path.dirname(args.msbayes_config_path),
            'pyco-nexus-files',
        )
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
    msbayes_cfg = pycoevolity.parsing.MsbayesConfig(
        args.msbayes_config_path,
        args.recode_ambig_states_as_missing,
    )
    sys.stdout.write('---\n')
    msbayes_cfg.write_ecoevolity_model_settings(out_stream = sys.stdout)
    sys.stdout.write("\ncomparisons:\n")
    for comparison_label, data in pycoevolity.parsing.Loci.iter_from_msbayes_config(
        msbayes_cfg,
        remove_triallelic_sites = args.remove_triallelic_sites,
        convert_to_binary = args.convert_to_binary,
    ):
        sys.stderr.write("Command: {0}\n".format(" ".join(argv)))
        sys.stderr.write("\tNumber of taxa:  {0}\n".format(data.number_of_taxa))
        sys.stderr.write("\tNumber of loci:  {0}\n".format(data.number_of_loci))
        sys.stderr.write("\tNumber of sites: {0}\n".format(data.number_of_sites))
        sys.stderr.write("\tNumber of triallelic sites found:   {0}\n".format(
                data.number_of_triallelic_sites_found))
        sys.stderr.write("\tNumber of triallelic sites removed: {0}\n".format(
                data.number_of_triallelic_sites_removed))

        out_name = f"{comparison_label}.nex"
        comp_path = f'{args.comparison_prefix}{out_name}'
        if not args.yaml_config_only:
            out_path = os.path.join(
                args.output_dir,
                out_name,
            )
            if os.path.exists(out_path):
                raise Exception("The path {0!r} already exists. "
                    "Please designate a different directory.".format(out_path))
            with open(out_path, "w") as out_stream:
                write_kwargs = {"stream": out_stream}
                if write_method == "write_nexus":
                    if args.charsets:
                        write_kwargs["include_charset_block"] = True
                getattr(data, write_method)(**write_kwargs)
        sys.stdout.write("- comparison:\n")
        sys.stdout.write(f"    path: \"{comp_path}\"\n")

def main_nexus(argv = sys.argv):
    main(argv = argv, write_method = "write_nexus")

def main_phylip(argv = sys.argv):
    main(argv = argv, write_method = "write_phylip")

if __name__ == "__main__":
    main()
