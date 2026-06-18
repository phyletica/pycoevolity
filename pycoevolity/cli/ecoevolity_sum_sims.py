#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import seaborn as sns

import pycoevolity


def parse_cli_args():
    parser = argparse.ArgumentParser(
        formatter_class = pycoevolity.argparse_utils.SmartDefaultsHelpFormatter,
    )

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
            'The number of parallel processes to use to parse results of '
            'analyses of simulated data sets.'
        ),
    )
    parser.add_argument(
        '-c', '--cred-interval-percent',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_percent_int,
        default = 95,
        help = (
            'The percent (as an integer) credible intervals to use '
            'for parameters.'
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
    parser.add_argument(
        '-d', '--model-distance',
        type = str,
        choices = ['edit', 'vi', 'ari'],
        default = 'edit',
        help = (
            'RAW|'
            'Choose one of three measures for calculating distances between\n'
            'the true event model and posterior samples of event models:\n'
            '  edit --- Edit distance; the minimal single-comparison transitions\n'
            '    to change one event model into another. Smaller values mean\n'
            '    samples are closer to the true model.\n'
            '  vi   --- Variation of Information distance. Smaller values mean\n'
            '    samples are closer to the true model. See\n'
            '    https://doi.org/10.1007/978-3-540-45167-9_14\n'
            '  ari  --- Adjusted Rand Index. Larger values mean samples are\n'
            '    closer to the true value. 1 is a perfect match, 0 is no better\n'
            '    than random, and negative values are worse than random. See\n'
            '    https://doi.org/10.1007/BF01908075\n'
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
        interval_percent = args.cred_interval_percent,
        model_distance_stat = args.model_distance,
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
    sys.stdout.write(
        f"Summary table of results written to '{summary_path}'\n"
    )

if __name__ == "__main__":
    main_cli()
