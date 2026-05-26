#!/usr/bin/env python

import os
import sys
import random
import argparse
import tempfile
import glob
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import pycoevolity
import pycoevolity.ecoevolity_config as eco_config
from .ecoevolity_viz_sims import write_existing_path_warning


def plot_cdf_comparison(
    prior_settings,
    posterior_samples,
    numpy_rng,
):
    if eco_config.distribution_is_fixed(prior_settings):
        prior_distribution = eco_config.get_fixed_distribution(prior_settings)
    else:
        prior_distribution = eco_config.sample_distribution(
            numpy_rng = numpy_rng,
            prior_settings = prior_settings,
            n = 100000,
        )
    fig = matplotlib.figure.Figure()
    gs = fig.add_gridspec(
        nrows = 1, ncols = 1,
        wspace = 0.0,
        hspace = 0.0,
    )
    ax = fig.add_subplot(gs[0, 0])
    empirical_line, model_line = pycoevolity.plotting.ax_compare_samples_to_cdf(
        ax = ax,
        samples = posterior_samples,
        prob_dist = prior_distribution,
        include_ks_test = True,
    )
    fig.tight_layout()
    return fig, ax, empirical_line, model_line

def plot_qq(
    prior_settings,
    posterior_samples,
    numpy_rng,
):
    if eco_config.distribution_is_fixed(prior_settings):
        prior_distribution = eco_config.get_fixed_distribution(prior_settings)
    else:
        prior_distribution = eco_config.sample_distribution(
            numpy_rng = numpy_rng,
            prior_settings = prior_settings,
            n = 100000,
        )
    fig = matplotlib.figure.Figure()
    gs = fig.add_gridspec(nrows = 1, ncols = 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = fig.add_subplot(gs[0, 0])
    qline = pycoevolity.plotting.ax_qq(
        ax = ax,
        samples = posterior_samples,
        prob_dist = prior_distribution,
    )
    fig.tight_layout()
    return fig, ax, qline

def process_parameter(
    parameter_name,
    parameter_settings,
    posterior_samples,
    output_prefix,
    numpy_rng,
    label = None,
    force = False,
    plot_ext = "svg",
):
    if not eco_config.parameter_is_estimated(parameter_settings):
        is_valid, expected_val, val = eco_config.fixed_param_values_are_valid(
            parameter_settings,
            posterior_samples,
        )
        if not is_valid:
            sys.stderr.write(
f"""
WARNING: According to the config file, '{parameter_name}' should be
fixed at {expected_val}, but value {val} was sampled.
Skipping plotting of '{parameter_name}'.
"""
            )
    else:
        plot_path = f"{output_prefix}prior-cdf-comparison-{parameter_name}.{plot_ext}"
        if (not force) and os.path.exists(plot_path):
            write_existing_path_warning(plot_path, sys.stderr)
            return
        if label is None:
            label = parameter_name
        fig, ax, eline, mline = plot_cdf_comparison(
            prior_settings = parameter_settings["prior"],
            posterior_samples = posterior_samples,
            numpy_rng = numpy_rng,
        )
        ax.set(xlabel = f"{label}")
        fig.savefig(plot_path, bbox_inches = "tight")
        plt.close(fig)
    
        plot_path = f"{output_prefix}prior-qq-plot-{parameter_name}.{plot_ext}"
        fig, ax, qline = plot_qq(
            prior_settings = parameter_settings["prior"],
            posterior_samples = posterior_samples,
            numpy_rng = numpy_rng,
        )
        ax.set(title = f"{label}")
        fig.savefig(plot_path, bbox_inches = "tight")
        plt.close(fig)

def process_event_model_prior(
    settings,
    nevent_samples,
    number_of_comparisons,
    output_prefix,
    numpy_rng,
    force = False,
    plot_ext = "svg",
):
    assert len(settings) == 1
    model_prior_name = list(settings.keys())[0]
    if model_prior_name == "fixed":
        return
    plot_path = f"{output_prefix}prior-cmf-comparison-nevents.{plot_ext}"
    if (not force) and os.path.exists(plot_path):
        write_existing_path_warning(plot_path, sys.stderr)
        return
    elif model_prior_name == "pitman_yor_process":
        model_prior_parameters = settings[model_prior_name][
                "parameters"]
        prior_nevent_samples = eco_config.sample_hyper_pitman_yor_distribution(
            numpy_rng = numpy_rng,
            prior_parameters = model_prior_parameters,
            number_of_elements = number_of_comparisons,
            n = 100000,
        )
    elif model_prior_name == "dirichlet_process":
        model_prior_parameters = settings[model_prior_name][
                "parameters"]
        prior_nevent_samples = eco_config.sample_hyper_dirichlet_distribution(
            numpy_rng = numpy_rng,
            prior_parameters = model_prior_parameters,
            number_of_elements = number_of_comparisons,
            n = 100000,
        )
    else:
        sys.stderr.write(
f"""
WARNING: Plotting of the prior sampling of the number of events for
event_model_prior '{model_prior_name}' is not currently supported.
Skipping plotting of the number of events.
"""
        )
        return
    fig = matplotlib.figure.Figure()
    gs = fig.add_gridspec(nrows = 1, ncols = 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = fig.add_subplot(gs[0, 0])
    emp_line, model_line = pycoevolity.plotting.ax_compare_nevents_samples_to_cdf(
        ax = ax,
        samples = nevent_samples,
        prior_samples = prior_nevent_samples,
        number_of_comparisons = number_of_comparisons,
        include_ks_test = True,
    )
    fig.tight_layout()
    fig.savefig(plot_path, bbox_inches = "tight")
    plt.close(fig)

def parse_cli_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        'config_path',
        metavar = 'ECOEVOLITY-CONFIG-PATH',
        type = pycoevolity.argparse_utils.arg_is_file,
        help = (
            'Path to ecoevolity configuration file to use to sample '
            'from the prior.'
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
        '-r', '--number-of-runs',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 10,
        help = (
            'The number of independent runs for sampling from the '
            'prior. This will determine the number of samples to be '
            'collected. The length of the chain and samping frequency '
            'are defined in the configuration file. The total number '
            'samples will equal (ignoring burn-in):  '
            'number of runs x ((chain_length / sample_frequency) - burnin).'
        ),
    )
    parser.add_argument(
        '-p', '--number-of-procs',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 4,
        help = (
            'The number of parallel processes to use to run the '
            '\'-r/--number-of-runs\' ecoevolity runs that sample from '
            'the prior. The default is the smaller of 4 or the number '
            'of runs.'
        ),
    )
    parser.add_argument(
        '-b', '--burnin',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_nonnegative_int,
        default = 0,
        help = (
            'The number of MCMC samples to remove from the beginning '
            'of each log file as burn in.'
        ),
    )
    parser.add_argument(
        '-s', '--step',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_positive_int,
        default = 1,
        help = (
            'The step value to use for thinning MCMC samples. '
            'For example, if set to 3, every third MCMC sample will be '
            'retained (other will be ignored).'
        ),
    )
    parser.add_argument(
        '--sparse-pop-size-plotting',
        action = 'store_true',
        help = ('Only plot pop sizes for first and last comparisons.'),
    )
    parser.add_argument(
        '-o', '--output-dir',
        action = 'store',
        default = os.curdir,
        type = pycoevolity.argparse_utils.arg_is_dir_or_new_dir,
        help = ('The directory in which to put all output files.'),
    )
    parser.add_argument(
        '--prefix',
        type = str,
        help = (
            'Prefix to add to the file name of every output file.'
        ),
    )
    parser.add_argument(
        '--plot-ext',
        type = str,
        default = 'svg',
        help = (
            'The file extension (and format) to use for output plotting files. '
            'Examples: '
            '\'--plot-ext svg\' (default), '
            '\'--plot-ext pdf\', '
            '\'--plot-ext png\', '
            '\'--plot-ext jpg\', etc. '
            'Any file formats supported by matplotlib should work.'
        ),
    )
    parser.add_argument(
        '--force',
        action = 'store_true',
        help = (
            'Overwrite output files if they already exist.'
        ),
    )
    parser.add_argument(
        '--seed',
        action = 'store',
        type = pycoevolity.argparse_utils.arg_is_positive_int,
        help = ('Seed for random number generator.'),
    )
    parser.add_argument(
        '--context',
        type = str,
        default = 'notebook',
        help = (
            'Value of \'context\' arguement passed to seaborn.sea_theme. '
            'See '
            'https://seaborn.pydata.org/generated/seaborn.set_theme.html '
            'for more info.'
        ),
    )
    parser.add_argument(
        '--style',
        type = str,
        default = 'ticks',
        help = (
            'Value of \'style\' arguement passed to seaborn.sea_theme. '
            'See '
            'https://seaborn.pydata.org/generated/seaborn.set_theme.html '
            'for more info.'
        ),
    )
    parser.add_argument(
        '--palette',
        type = str,
        default = 'colorblind',
        help = (
            'Value of \'palette\' arguement passed to seaborn.sea_theme. '
            'See '
            'https://seaborn.pydata.org/generated/seaborn.set_theme.html '
            'for more info.'
        ),
    )
    parser.add_argument(
        '--font',
        type = str,
        default = 'sans-serif',
        help = (
            'Value of \'font\' arguement passed to seaborn.sea_theme. '
            'See '
            'https://seaborn.pydata.org/generated/seaborn.set_theme.html '
            'for more info.'
        ),
    )
    parser.add_argument(
        '--font-scale',
        type = pycoevolity.argparse_utils.arg_is_positive_float,
        default = 1.0,
        help = (
            'Value of \'font-scale\' arguement passed to seaborn.sea_theme. '
            'See '
            'https://seaborn.pydata.org/generated/seaborn.set_theme.html '
            'for more info.'
        ),
    )
    return parser.parse_args()

def check_for_existing_outputs(output_prefix, file_ext):
    output_pattern = f"{output_prefix}prior-*.{file_ext}"
    output_files = glob.glob(output_pattern)
    if output_files:
        msg = (
            "ERROR: The following existing files would be overwritten.\n"
            "If you wish to overwrite them, please use the \'--force\' "
            "option.\n\t{0}\n".format(
                "\n\t".join(output_files)
            )
        )
        sys.stderr.write(msg)
        sys.exit(1)

def main_cli():
    args = parse_cli_args()

    sns.set_theme(
        context = args.context,
        style = args.style,
        palette = args.palette,
        font = args.font,
        font_scale = args.font_scale,
    )

    output_dir = pycoevolity.argparse_utils.process_output_dir_arg(args.output_dir)
    prefix = ''
    if args.prefix:
        prefix = args.prefix
    output_prefix = os.path.join(
        output_dir,
        prefix,
    )

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
    np_rng = pycoevolity.rng_utils.get_numpy_rng(seed)

    config_name = os.path.splitext(os.path.basename(args.config_path))[0]
    output_prefix = f"{output_prefix}{config_name}-"

    if not args.force:
        check_for_existing_outputs(output_prefix, args.plot_ext)

    seeds = pycoevolity.rng_utils.get_safe_seeds(rng, n = args.number_of_runs)

    tmp_prefix = "ecoevolity-prior-sampling-"
    with tempfile.TemporaryDirectory(prefix = tmp_prefix) as temp_dir:
        sys.stdout.write(
            "Collecting ecoevolity MCMC samples from the prior distribution...\n"
        )
        state_log_paths = pycoevolity.ecoevolity.collect_prior_samples(
            seeds = seeds,
            config_path = args.config_path,
            output_dir = temp_dir,
            number_of_procs = args.number_of_procs,
            working_dir = None,
            eco_exe_dir = eco_exe_dir,
            timeout = 600,
            max_num_attempts = 2,
        )
        sys.stdout.write(
            "Parsing MCMC samples...\n"
        )
        posterior_sample = pycoevolity.posterior.PosteriorSample(
            paths = state_log_paths,
            burnin = args.burnin,
            step = args.step,
            include_time_in_coal_units = False,
        )
        sys.stdout.write(
            f"Total number of MCMC samples retained: "
            f"{posterior_sample.number_of_samples}\n"
        )

    config = eco_config.get_yaml_config(args.config_path)
    number_of_comparisons = len(config.get("comparisons"))

    model_prior_settings = config.get("event_model_prior")
    assert len(model_prior_settings) == 1

    sys.stdout.write(
        "Plotting the number of events...\n"
    )
    process_event_model_prior(
        settings = model_prior_settings,
        nevent_samples = posterior_sample.parameter_samples["number_of_events"],
        number_of_comparisons = number_of_comparisons,
        output_prefix = output_prefix,
        numpy_rng = np_rng,
        force = args.force,
        plot_ext = args.plot_ext,
    )

    model_prior_name = list(model_prior_settings.keys())[0]
    if not model_prior_name == "fixed":
        model_prior_parameters = model_prior_settings[model_prior_name][
                "parameters"]
        for parameter_name, parameter_settings in model_prior_parameters.items():
            values = posterior_sample.parameter_samples[parameter_name]
            label = parameter_name
            if model_prior_name == "dirichlet_process":
                label = f"Dirichlet process {parameter_name}"
            elif model_prior_name == "pitman_yor_process":
                label = f"Pitman Yor process {parameter_name}"
            process_parameter(
                parameter_name = parameter_name,
                parameter_settings = parameter_settings,
                posterior_samples = values,
                output_prefix = output_prefix,
                numpy_rng = np_rng,
                label = label,
                force = args.force,
                plot_ext = args.plot_ext,
            )

    event_time_prior = config["event_time_prior"]
    assert len(event_time_prior) == 1
    event_time_prior_name = list(event_time_prior.keys())[0]

    sys.stdout.write(
        "Plotting the event time prior parameters (if any)...\n"
    )
    event_time_prior_parameters = event_time_prior[event_time_prior_name]
    num_prior_parameters = len(event_time_prior_parameters)
    for parameter, settings in event_time_prior_parameters.items():
        if eco_config.parameter_is_estimated(settings):
            # TODO: This is a brittle hack to guess the header for this
            # parameter in the ecoevolity state log output. The risk of
            # guessing wrong is low, because if we get it wrong the prior
            # sampling will look very wrong (i.e., guessing wrong will never
            # mask a real problem with MCMC sampling).
            param_key_idx = 0
            if (
                (num_prior_parameters > 1) and
                (parameter in ("max", "scale", "standard_deviation"))
            ):
                param_key_idx = 1
            param_key = f"time_prior_parameter_{param_key_idx}"
            values = posterior_sample.parameter_samples[param_key]
            label = f"Time prior {parameter}"
            process_parameter(
                parameter_name = parameter,
                parameter_settings = settings,
                posterior_samples = values,
                output_prefix = output_prefix,
                numpy_rng = np_rng,
                label = label,
                force = args.force,
                plot_ext = args.plot_ext,
            )

    event_time_settings = {
        "estimate" : True,
        "prior" : event_time_prior,
    }

    sys.stdout.write(
        "Plotting the event times...\n"
    )
    height_keys = list(posterior_sample.get_height_keys())
    # Vet div times for first and last comparison
    height_keys_to_plot = [height_keys[0], height_keys[-1]]
    for height_key in height_keys_to_plot:
        values = posterior_sample.parameter_samples[height_key]
        process_parameter(
            parameter_name = height_key,
            parameter_settings = event_time_settings,
            posterior_samples = values,
            output_prefix = output_prefix,
            numpy_rng = np_rng,
            force = args.force,
            plot_ext = args.plot_ext,
        )

    ##########################################################################
    # Handle population sizes parameters
    ##########################################################################

    sys.stdout.write(
        "Plotting the population sizes...\n"
    )

    global_comp_settings = config.get("global_comparison_settings", {})

    comps_to_plot = config["comparisons"]
    if args.sparse_pop_size_plotting:
        comps_to_plot = [
            config["comparisons"][0],
            config["comparisons"][-1],
        ]

    # Handle leaf population size parameters
    for comp_idx, comp in enumerate(comps_to_plot):
        comp_settings = global_comp_settings
        comp_settings.update(comp["comparison"])
        # It's only worth processing the parameter if we have prior settings
        # for it. We don't want to try and guess what the defaults are for
        # ecoevolity
        if ("parameters" in comp_settings) and ("population_size" in comp_settings["parameters"]):
            # Only process first leaf population
            leaf_label = posterior_sample.tip_labels[comp_idx][0]
            leaf_pop_size_key = f"pop_size_{leaf_label}"
            values = posterior_sample.parameter_samples[leaf_pop_size_key]
            process_parameter(
                parameter_name = leaf_pop_size_key,
                parameter_settings = comp_settings["parameters"]["population_size"],
                posterior_samples = values,
                output_prefix = output_prefix,
                numpy_rng = np_rng,
                force = args.force,
                plot_ext = args.plot_ext,
            )
    # Handle root population size parameters
    for comp_idx, comp in enumerate(comps_to_plot):
        comp_settings = global_comp_settings
        comp_settings.update(comp["comparison"])
        comp_equal_pop_sizes = comp_settings["equal_population_sizes"]
        root_pop_size_key = f"pop_size_root_{posterior_sample.height_labels[comp_idx]}"
        if comp_equal_pop_sizes:
            # No need to plot this root pop size, because it's identical to the
            # leaf pop size plotted above. However, we will verify that the
            # values are identical
            leaf_pop_size_key = f"pop_size_{posterior_sample.height_labels[comp_idx]}"
            leaf_sizes = posterior_sample.parameter_samples[leaf_pop_size_key]
            root_sizes = posterior_sample.parameter_samples[root_pop_size_key]
            assert len(leaf_sizes) == len(root_sizes)
            for sample_idx in range(len(leaf_sizes)):
                if not pycoevolity.math_utils.diff_almost_zero(
                    leaf_sizes[sample_idx],
                    root_sizes[sample_idx],
                ):
                    msg = (
f"""
Comparison {comp_idx + 1} is set to have equal population sizes, but values
from Sample {sample_idx + 1} are not equal:
    Root size: {root_sizes[sample_idx]}
    Leaf size: {leaf_sizes[sample_idx]}
"""
                    )
                    raise Exception(msg)
        else:
            # Process unconstrained root pop size if we have prior settings for it
            if ("parameters" in comp_settings) and (
                "root_relative_population_size" in comp_settings["parameters"]):
                leaf_pop_size_keys = [
                    f"pop_size_{l}" for l in posterior_sample.tip_labels[comp_idx]
                ]
                values = []
                # Convert root pop sizes to relative sizes to be on the same scale
                # as the prior
                for sample_idx in range(posterior_sample.number_of_samples):
                    root_size = posterior_sample.parameter_samples[root_pop_size_key][sample_idx]
                    leaf_sizes = [
                        posterior_sample.parameter_samples[k][sample_idx] for k in leaf_pop_size_keys
                    ]
                    assert (len(leaf_sizes) == 1) or (len(leaf_sizes) == 2)
                    mean_leaf_size = sum(leaf_sizes) / len(leaf_sizes)
                    rel_root_size = root_size / mean_leaf_size
                    values.append(rel_root_size)
                process_parameter(
                    parameter_name = root_pop_size_key,
                    parameter_settings = comp_settings["parameters"]["root_relative_population_size"],
                    posterior_samples = values,
                    output_prefix = output_prefix,
                    numpy_rng = np_rng,
                    force = args.force,
                    plot_ext = args.plot_ext,
                )

if __name__ == "__main__":
    main_cli()
