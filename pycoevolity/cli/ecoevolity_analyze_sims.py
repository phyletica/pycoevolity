#!/usr/bin/env python

import os
import sys
import argparse
import signal

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
        '-e', '--ecoevolity-dir',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_dir,
        help = (
            'The directory in which ecoevolity\'s programs are '
            'installed. By default, ecoevolity\'s programs will be ' 
            'called without a path (i.e., the directory in which '
            'they are installed need to be in your environment\'s '
            'PATH variable.'
        ),
    )
    parser.add_argument(
        '-p', '--number-of-procs',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 4,
        help = (
            'The number of parallel processes to use to analyze '
            'simulations with ecoevolity.'
        ),
    )
    parser.add_argument(
        '-t', '--chain-timeout',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 86400,
        help = (
            'The timeout (in seconds) for each ecoevolity MCMC chain. '
            'If a chain runs longer than this the subprocess will stop the '
            'chain and raise an error.'
        ),
    )
    parser.add_argument(
        '-a', '--chain-attempts',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 3,
        help = (
            'The max number of times to try running each ecoevolity MCMC '
            'chain.'
        ),
    )

    args = parser.parse_args()
    return args

def main_cli():
    args = parse_cli_args()

    eco_exe_dir = os.path.abspath(
        pycoevolity.ecoevolity.get_ecoevolity_dir(
            dir_to_check = args.ecoevolity_dir,
        )
    )

    sys.stdout.write(
        f"Using ecoevolity programs found in '{eco_exe_dir}'\n"
    )

    json_path = os.path.abspath(args.json_path)
    json_dir = os.path.dirname(json_path)
    results = pycoevolity.fileio.load_json(json_path)
    rel_sim_files_dir = results["simulation_files_dir"]

    # If the process is killed externally, let's make sure the analyses that
    # have finished get recorded in the json file
    def handle_exit(signum, frame):
        pycoevolity.fileio.write_json(results, json_path, indent = 4)
    signal.signal(signal.SIGTERM, handle_exit)
    signal.signal(signal.SIGINT, handle_exit)

    sys.stdout.write(
        f"Starting ecoevolity analyses of simulated datasets...\n"
    )
    num_analyses = 0
    try:
        num_analyses = pycoevolity.ecoevolity.run_analyses_on_sims(
            sim_data = results["simulations"],
            working_dir = json_dir,
            eco_exe_dir = eco_exe_dir,
            number_of_procs = args.number_of_procs,
            relax_constant_sites = False,
            relax_missing_sites = False,
            relax_triallelic_sites = False,
            timeout = args.chain_timeout,
            max_num_attempts = args.chain_attempts,
            output_dir = rel_sim_files_dir,
        )
    except Exception as e:
        # Write results so any completed analyses won't need to be re-run
        pycoevolity.fileio.write_json(results, json_path, indent = 4)
        raise e

    if num_analyses == 0:
        sys.stdout.write(
            "All ecoevolity analyses of simulated data sets were already "
            "complete.\n"
        )
    else:
        # Write results now, so we don't have to repeat ecoevolity analyses if
        # something goes wrong with sumcoevolity below
        pycoevolity.fileio.write_json(results, json_path, indent = 4)
    sys.stdout.write(
        f"Running sumcoevolity on ecoevolity results...\n"
    )
    num_analyses = 0
    try:
        num_analyses = pycoevolity.ecoevolity.add_sumcoevolity_to_results(
            sim_data = results["simulations"],
            working_dir = json_dir,
            eco_exe_dir = eco_exe_dir,
            output_dir = rel_sim_files_dir,
            num_prior_draws = results["number_of_prior_draws"],
            burnin = results["burnin"],
            number_of_procs = args.number_of_procs,
        )
    except Exception as e:
        # Write results so any completed analyses won't need to be re-run
        pycoevolity.fileio.write_json(results, json_path, indent = 4)
        raise e
    if num_analyses == 0:
        sys.stdout.write(
            "All sumcoevolity analyses of ecoevolity results were already "
            "complete.\n"
        )
    else:
        pycoevolity.fileio.write_json(results, json_path, indent = 4)

if __name__ == "__main__":
    main_cli()
