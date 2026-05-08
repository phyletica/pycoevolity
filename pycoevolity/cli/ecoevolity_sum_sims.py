#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import seaborn as sns

import pycoevolity


def parse_cli_args():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'json_path',
        metavar = 'SIM-DATA-JSON-FILE',
        type = pycoevolity.argparse_utils.arg_is_file,
        help = (
            'Path to JSON-formatted simulation data file.'
        ),
    )
    parser.add_argument(
        '-p', '--number-of-procs',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 4,
        help = (
            'The number of processors to use to parse results of analyses '
            'of simulated data sets.'
        ),
    )
    parser.add_argument(
        '-l', '--config-label-file',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_file,
        help = (
            'Path to YAML-formatted file that maps config file names to labels '
            'used for plotting.'
        ),
    )
    parser.add_argument(
        '--include-time-in-coal-units',
        action = 'store_true',
        help = (
            'Include results of event times converted to coalescent units.'
        ),
    )

    args = parser.parse_args()
    return args

def main_cli():
    args = parse_cli_args()
    config_labels = None
    if args.config_label_file:
        config_labels = pycoevolity.parsing.parse_config_label_yaml(
            args.config_label_file)
    df = pycoevolity.ecoevolity.parse_sim_results(
        results_path = args.json_path,
        config_labels = config_labels,
        include_time_in_coal_units = args.include_time_in_coal_units,
        number_of_procs = args.number_of_procs,
    )
    summary_path = os.path.join(
        os.path.dirname(args.json_path),
        "results-summary.tsv.gz",
    )
    df.to_csv(
        summary_path,
        sep = "\t",
        compression = "gzip",
        index = False,
    )

if __name__ == "__main__":
    main_cli()
