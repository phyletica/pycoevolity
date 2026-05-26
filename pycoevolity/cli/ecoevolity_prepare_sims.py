#!/usr/bin/env python

import os
import sys
import random
import argparse

import pycoevolity


def parse_cli_args():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'config_paths',
        metavar = 'ECOEVOLITY-CONFIG-PATH',
        type = pycoevolity.argparse_utils.arg_is_file,
        nargs = "*",
        help = (
            'Paths to ecoevolity configuration files to use to analyze each '
            'simulated data set.'
        ),
    )
    parser.add_argument(
        '-c', '--sim-config',
        metavar = 'ECOEVOLITY-CONFIG-PATH',
        dest = 'sim_configs',
        action = 'append',
        required = False,
        type = pycoevolity.argparse_utils.arg_is_file,
        help = (
            'Path to the ecoevolity configuration file to use to simulate '
            'data sets from the prior. This option can be used multiple times '
            'if you want to generate simulations under multiple configs.'
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
        '-s', '--number-of-sims',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 100,
        help = (
            'The number of data sets to simulate with simcoevolity.'
        ),
    )
    parser.add_argument(
        '-r', '--number-of-chains',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 2,
        help = (
            'The number of independent ecoevolity MCMC chains '
            'to run on each siumlated data set.'
        ),
    )
    parser.add_argument(
        '-p', '--number-of-procs',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 4,
        help = (
            'The number of parallel processes to use to simulate data sets '
            'with simcoevolity.'
            'NOTE: changing the number of processes will change the results '
            'for a given random seed. '
            'Thus to reproduce identical results you need to specify the '
            'same seed AND the same number of processes using the '
            '\'-p\'/\'--number-of-procs\' argument.'
        ),
    )
    parser.add_argument(
        '-b', '--burnin',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 101,
        help = (
            'The number of samples to remove from the beginning of '
            'each log file as burn in.'
        ),
    )
    parser.add_argument(
        '--number-of-prior-draws',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 100000,
        help = (
            'The number of prior samples to use when running sumcoevolity.'
        ),
    )
    parser.add_argument(
        '-o', '--output-dir',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_dir_or_new_dir,
        help = (
            'The directory in which to put all output files.'
        ),
    )
    parser.add_argument(
        '--seed',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_positive_int,
        help = (
            'Seed for random number generator. '
            'NOTE: To reproduce identical results, you need to specify the '
            'same seed AND the same number of processes using the '
            '\'-p\'/\'--number-of-procs\' argument.'
        ),
    )
    parser.add_argument(
        '--append-to',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_file,
        help = (
            'Path to results file to which to append more simulations. '
            'When this option is used, the only other options that are used '
            'are: '
            '--seed, '
            '-e | --ecoevolity-dir, '
            '-s | --number-of-sims, '
            'and '
            '-p | --number-of-procs, '
            'All other arguments will be ignored.'
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

    rng = random.Random()
    seed = args.seed
    if not seed:
        seed = pycoevolity.rng_utils.get_safe_seed()
    rng.seed(seed)

    results = None
    json_path = None
    json_dir = None
    output_dir = None
    sim_files_dir = None
    if args.append_to:
        json_path = args.append_to
        json_dir = os.path.abspath(os.path.dirname(json_path))
        results = pycoevolity.fileio.load_json(json_path)
        if seed in results["seeds"]:
            raise Exception(
                f"Seed {seed} was already used; please use a different seed."
            )
        results["seeds"].append(seed)
        args.sim_configs = [
            os.path.abspath(
                os.path.join(json_dir, p)
            ) for p in results["simulation_configs"]
        ]
        args.config_paths = [
            os.path.abspath(
                os.path.join(json_dir, p)
            ) for p in results["inference_configs"]
        ]
        args.number_of_chains = results["number_of_chains"]
        args.number_of_prior_draws = results["number_of_prior_draws"]
        args.burnin = results["burnin"]
        sim_files_dir = os.path.abspath(
            os.path.join(
                json_dir,
                results["simulation_files_dir"],
            )
        )
        output_dir = json_dir
        if not os.path.isdir(sim_files_dir):
            raise Exception(
                f"Unexpected results directory determined from json results "
                f"file: {sim_files_dir}\n"
                "Please don't move a 'simulation-data.json' file before appending "
                "simulations to it."
            )
    else:
        if (not args.sim_configs) or (not args.config_paths):
            raise Exception(
                "Simulation and inference configs are required when not "
                "appending to previous results."
            )
        args.sim_configs = [os.path.abspath(p) for p in args.sim_configs]
        args.config_paths = [os.path.abspath(p) for p in args.config_paths]
        output_dir = os.path.abspath(
            pycoevolity.argparse_utils.process_output_dir_arg(args.output_dir))
        sim_files_dir = pycoevolity.argparse_utils.process_output_dir_arg(
            os.path.join(output_dir, "simulation-files")
        )
        json_path = os.path.join(output_dir, "simulation-data.json")
        if os.path.exists(json_path):
            raise Exception(
f"""
Simulation data file already exists: \'{json_path}\'
If you want to append more simulations, please use the \'--append-to\' argument.
To learn more about the \'--append-to\' argument you can use the \'-h\' or
\'--help\' flags to see the help menu.
"""
            )
        json_dir = output_dir
        results = { "seeds" : [seed] }
        results["simulation_configs"] = [
            os.path.relpath(p, json_dir) for p in args.sim_configs
        ]
        results["inference_configs"] = [
            os.path.relpath(p, json_dir) for p in args.config_paths
        ]
        results["number_of_chains"] = args.number_of_chains
        results["number_of_prior_draws"] = args.number_of_prior_draws
        results["burnin"] = args.burnin
        results["simulation_files_dir"] = os.path.relpath(sim_files_dir, json_dir)

    if not pycoevolity.fileio.file_names_are_unique(args.sim_configs):
        raise Exception(
            "Simulation config file names are not unique"
        )
    if not pycoevolity.fileio.file_names_are_unique(args.config_paths):
        raise Exception(
            "Inference config file names are not unique"
        )
    if not pycoevolity.ecoevolity_config.configs_share_comparisons(
        yaml_config_paths = set(args.sim_configs + args.config_paths),
    ):
        raise Exception(
            "All configs do not share the same comparisons"
        )

    sim_configs = results["simulation_configs"]
    infer_configs = results["inference_configs"]
    rel_sim_files_dir = os.path.relpath(sim_files_dir, json_dir)

    sys.stdout.write(
        f"Generating simulated datasets...\n"
    )
    sim_results = pycoevolity.ecoevolity.prepare_simulations(
        rng = rng,
        sim_configs = sim_configs,
        infer_configs = infer_configs,
        output_dir = rel_sim_files_dir,
        working_dir = json_dir,
        eco_exe_dir = eco_exe_dir,
        number_of_sims = args.number_of_sims,
        number_of_procs = args.number_of_procs,
        number_of_chains = args.number_of_chains,
        singleton_sample_prob = None,
        locus_size = None,
        max_one_variable_site_per_locus = False,
        charsets = False,
        relax_constant_sites = False,
        relax_missing_sites = False,
        relax_triallelic_sites = False,
        output_nexus = False,
    )
    if args.append_to:
        append_results(results, sim_results)
    else:
        results["simulations"] = sim_results
    pycoevolity.fileio.write_json(results, json_path, indent = 4)
    sys.stdout.write(
        f"Simulation data written to '{json_path}'\n"
        f"Simulation files output in '{sim_files_dir}'\n"
    )

if __name__ == "__main__":
    main_cli()
